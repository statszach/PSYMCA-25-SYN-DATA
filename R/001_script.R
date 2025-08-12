#### Load packages

library(tidyverse)
library(synthpop)
library(MplusAutomation)
library(fs)
library(here)
library(gtsummary)

#### Read data

# Zach's file path
HCAP2016 <- readRDS(fs::path("C:", "Users", "zkunicki", 
                               "Brown Dropbox", "Zachary Kunicki",
                               "PsyMCA2025_Sharing", "Data",
                               "HRS-HCAP-Classification",
                               "HCAPHRS.RDS"))

#### Pull vars of interest

demos_obs <- HCAP2016 %>% 
  dplyr::filter(age65up == 1) %>% 
  dplyr::select(rage, female, black, hisp, SCHLYRS)

cog_obs <- HCAP2016 %>% 
  dplyr::filter(age65up == 1) %>% 
  dplyr::select(vdori, vdlfl1z,
                vdlfl2, vdlfl3,
                vdwdimmz, vdwddelz,
                vdexf7z, vdsevens,
                vdcount)

#### Generate synthetic datasets

demos_synth <- synthpop::syn(demos_obs, seed = '0904', minnumlevels = 2)

cog_synth <- synthpop::syn(cog_obs, seed = '0904', minnumlevels = 5)

demos_syn <- demos_synth$syn %>% 
  dplyr::mutate(group = "Synthesized")

cog_syn <- cog_synth$syn %>% 
  dplyr::mutate(group = "Synthesized")

#### Make long-form datasets via smushing

demos_merged <- demos_obs %>% 
  dplyr::mutate(group = "Observed") %>% 
  bind_rows(demos_syn)

cog_merged <- cog_obs %>% 
  dplyr::mutate(group = "Observed") %>% 
  bind_rows(cog_syn)

#### Check for differences using gtsummary

demos_merged %>% 
  gtsummary::tbl_summary(by = group,
                         statistic = list(all_continuous() ~ "{mean} ({sd})")) %>% 
  gtsummary::add_p(test = list(all_continuous() ~ "t.test")) 

cog_merged %>% 
  gtsummary::tbl_summary(by = group,
                         statistic = list(all_continuous() ~ "{mean} ({sd})")) %>% 
  gtsummary::add_p(test = list(all_continuous() ~ "t.test")) 