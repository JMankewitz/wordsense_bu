#take in a partition of data and run a model with bootstrapped confidence intervals

library("lme4")
library("lmerTest")
#library("optimx")
library("merTools")
library("boot")
library('tictoc')
library('here')
library('tidyverse')


i_am_finished = function(){
  command = paste0("osascript -e 'display notification \"R has finished the requested processing task.\" with title \"Processing complete!\"'")
  system(command)
}

message("Process Data")


full_sense_tags <- readRDS(here("data/processed_data/full_sense_tags.rds"))

run_model_for_speaker <- function(target_speaker_role, start_age, stop_age, bootstrap_num){
  maj_tags <- full_sense_tags %>%
    #drop other meanings
    filter(sense_name != "other_meanings") %>% 
    #restrict to type+sense combinations that were seen at least once - drop noise from unlikely tags
    group_by(target_child_name, type_sense) %>%
    mutate(n_type_sense = n()) %>% 
    filter(n_type_sense > 1) %>% ungroup() %>%
    #restrict to tags in target partition
    filter(speaker_role == target_speaker_role,
           start_age_in_months >= start_age,
           start_age_in_months < stop_age,
           !sense_name %in% c("wrong_pos", "other_meanings"))
  
  #floor start age
  maj_tags$downsampled_start_age_in_months = sapply(maj_tags$start_age_in_months, 
                                                    function(x){floor(x /6)*6})
  
  maj_tags$downsampled_start_age_in_months_adjusted = maj_tags$downsampled_start_age_in_months - start_age
  
  #collect the most frequent senses
  most_freq_senses <- aggregate(sense_name ~ type + downsampled_start_age_in_months + target_child_name, 
                                maj_tags, 
                                function(x){x_tab = table(x)
                                return(names(x_tab[order(x_tab)])[1])
                                }
  ) %>% rename("most_common_sense_name" = "sense_name")
  
  maj_tags <- maj_tags %>% merge(most_freq_senses) %>%
    #1 if the sense name is anything except the most common sense in that interval- represents polysemy!
    mutate(is_most_common_sense_name = as.numeric(sense_name == most_common_sense_name),
           is_not_most_common_sense_name = as.numeric(sense_name != most_common_sense_name))
  
  #run model
  message("Run Model")
  
  model = glmer(is_not_most_common_sense_name ~  downsampled_start_age_in_months_adjusted  + (downsampled_start_age_in_months_adjusted | target_child_name), 
                data = subset(maj_tags, !is.na(is_most_common_sense_name)), 
                family = binomial,
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=400000)))
  
  #bootstrap
  new_data <- maj_tags %>% distinct(downsampled_start_age_in_months_adjusted)
  
  message("Get CIs")
  
  model_boot <- bootMer(model, nsim=bootstrap_num, 
                        re.form=NA,  #.progress = "txt", 
                        parallel = "multicore",
                        ncpus = 2,
                        FUN = function(x)inv.logit(predict(x, newdata = new_data, re.form=NA)))
  
  saveRDS(model_boot, paste0(here("data/derived_data/"), target_speaker_role, "_model_full.rds"))
  
  get_ci <- function(x, quant){quantile(x, probs = quant)}
  
  low_ci <- apply(model_boot$t, 2, get_ci, quant = .025)
  high_ci <- apply(model_boot$t, 2, get_ci, quant = 0.975)
  ci_model_data <- cbind(maj_tags %>% distinct(grouped_age_interval, downsampled_start_age_in_months_adjusted), low_ci) %>% 
    mutate(speaker_role = target_speaker_role) %>% cbind(high_ci)
  #i_am_finished()
  saveRDS(ci_model_data, paste0(here("data/derived_data/"), target_speaker_role, "_model_CIs.rds"))
}

run_model_for_speaker(target_speaker_role = "Child",
                      start_age = 12, 
                      stop_age = 48,
                      bootstrap_num = 1000)

run_model_for_speaker(target_speaker_role = "Caregiver",
                      start_age = 12, 
                      stop_age = 48,
                      bootstrap_num = 1000)
