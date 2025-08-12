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

demo_and_cog <- HCAP2016 %>% 
  dplyr::filter(age65up == 1) %>% 
  dplyr::select(rage, female, black, hisp, SCHLYRS,
                vdori, vdlfl1z,
                vdlfl2, vdlfl3,
                vdwdimmz, vdwddelz,
                vdexf7z, vdsevens,
                vdcount)


#### Generate synthetic datasets

cog_and_demo_synth <- synthpop::syn(demo_and_cog, seed = '0904', minnumlevels = 5)

cog_and_demo_syn <- cog_and_demo_synth$syn %>%
  dplyr::mutate(group = "Synthesized")

#### Make long-form datasets via smushing

demo_and_cog_merged <- demo_and_cog %>% 
  dplyr::mutate(group = "Observed") %>% 
  dplyr::bind_rows(cog_and_demo_syn)

haven::write_dta(demo_and_cog_merged,
                 "demo_and_cog_merged.dta")

#### Check for differences using gtsummary

demo_and_cog_merged %>% 
  gtsummary::tbl_summary(by = group,
                         statistic = list(all_continuous() ~ "{mean} ({sd})")) %>% 
  gtsummary::add_p(test = list(all_continuous() ~ "t.test")) 

# No differences emerged -- but the absence of evidence is not evidence of absence.

#### Develop factor model of cognition in the observed data

singledomainmodel_obs <- mplusObject(
  MODEL = "
  gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
  gcp@1;
  
  recall BY vdwdimmz* (1);
  recall BY vdwddelz* (1);
  recall@1;
  
  gcp with recall@0;
  

  ",
  usevariables = colnames(cog_obs),
  VARIABLE = "CATEGORICAL = vdori vdlfl2 vdlfl3
  vdsevens vdcount;",
  OUTPUT = "TECH4 STDYX CINTERVAL MODINDICES(ALL);",
  ANALYSIS = "ESTIMATOR = WLSMV;
  PARAMETERIZATION = THETA;",
  rdata = cog_obs)

gcp_obs <- mplusModeler(singledomainmodel_obs, 
                    modelout = "gcpobs.inp", run = TRUE)

#### Run factor model in synthetic data

synmplus <- cog_syn %>% dplyr::select(-group)

singledomainmodel_syn <- mplusObject(
  MODEL = "
  gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
  gcp@1;
  
  recall BY vdwdimmz* (1);
  recall BY vdwddelz* (1);
  recall@1;
  
  gcp with recall@0;
  

  ",
  usevariables = colnames(synmplus),
  VARIABLE = "CATEGORICAL = vdori vdlfl2 vdlfl3
  vdsevens vdcount;",
  OUTPUT = "TECH4 STDYX CINTERVAL MODINDICES(ALL);",
  ANALYSIS = "ESTIMATOR = WLSMV;
  PARAMETERIZATION = THETA;",
  rdata = synmplus)

gcp_syn <- mplusModeler(singledomainmodel_syn, 
                        modelout = "gcpsyn.inp", run = TRUE)

#### Move to invariance models

configural_model <- mplusObject(
  MODEL = "
  
  gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
  gcp@1;
  
  recall BY vdwdimmz* (1);
  recall BY vdwddelz* (1);
  recall@1;
  
  gcp with recall@0;
  
  MODEL OBSERVED:
  
  gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
  gcp@1;
  
  recall BY vdwdimmz* (2);
  recall BY vdwddelz* (2);
  recall@1;
  
  gcp with recall@0;

  [ vdori$1*] ;
  [ vdori$2*] ;
  [ vdori$3*] ;
  [ vdori$4*] ;
  [ vdlfl2$1*] ;
  [ vdlfl2$2*] ;
  [ vdlfl3$1*] ;
  [ vdlfl3$2*] ;
  [ vdsevens$1*] ;
  [ vdsevens$2*] ;
  [ vdsevens$3*] ;
  [ vdsevens$4*] ;
  [ vdsevens$5*] ;
  [ vdcount$1*] ;
  
  [ vdlfl1z*] ;
  [ vdwdimmz*] ;
  [ vdwddelz* ] ;
  [ vdexf7z* ] ;
  
  MODEL SYNTHETIC:
  
  gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
  gcp@1;
  
  recall BY vdwdimmz* (3);
  recall BY vdwddelz* (3);
  recall@1;
  
  gcp with recall@0;
  
  [ vdori$1*] ;
  [ vdori$2*] ;
  [ vdori$3*] ;
  [ vdori$4*] ;
  [ vdlfl2$1*] ;
  [ vdlfl2$2*] ;
  [ vdlfl3$1*] ;
  [ vdlfl3$2*] ;
  [ vdsevens$1*] ;
  [ vdsevens$2*] ;
  [ vdsevens$3*] ;
  [ vdsevens$4*] ;
  [ vdsevens$5*] ;
  [ vdcount$1*] ;
  
  [ vdlfl1z*] ;
  [ vdwdimmz*] ;
  [ vdwddelz* ] ;
  [ vdexf7z* ] ;
  
  [gcp@0];
  [recall@0];
  vdori@1;
  vdlfl2@1;
  vdlfl3@1;
  vdsevens@1;
  vdcount@1;

  

  ",
  usevariables = colnames(cog_merged),
  VARIABLE = "CATEGORICAL = vdori vdlfl2 vdlfl3
  vdsevens vdcount;
  GROUP = group (1 = Observed 2 = Synthetic);",
  OUTPUT = "TECH4 STDYX SVALUES;",
  ANALYSIS = "ESTIMATOR = WLSMV;
  PARAMETERIZATION = THETA;",
  rdata = cog_merged)

configural <- mplusModeler(configural_model, 
                        modelout = "configural.inp", run = TRUE)

#### Metric model (constrain loadings and intercepts)

metric_model <- mplusObject(
  MODEL = "
  
  gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
  gcp@1;
  
  recall BY vdwdimmz* (1);
  recall BY vdwddelz* (1);
  recall@1;
  
  gcp with recall@0;
  
  MODEL OBSERVED:
  
  gcp BY vdori* (l1);
  gcp BY vdlfl1z (l2);
  gcp BY vdlfl2 (l3); 
  gcp BY vdlfl3 (l4);
  gcp BY vdwdimmz (l5);
  gcp BY vdwddelz (l6);
  gcp BY vdexf7z (l7);
  gcp BY vdsevens (l8);
  gcp BY vdcount (l9);
  gcp@1;
  
  recall BY vdwdimmz* (ml1);
  recall BY vdwddelz* (ml1);
  recall@1;
  
  gcp with recall@0;

  [ vdori$1*] (t1);
  [ vdori$2*] (t2);
  [ vdori$3*] (t3);
  [ vdori$4*] (t4);
  [ vdlfl2$1*] (t5);
  [ vdlfl2$2*] (t6);
  [ vdlfl3$1*] (t7);
  [ vdlfl3$2*] (t8);
  [ vdsevens$1*] (t9);
  [ vdsevens$2*] (t10);
  [ vdsevens$3*] (t11);
  [ vdsevens$4*] (t12);
  [ vdsevens$5*] (t13);
  [ vdcount$1*] (t14);
  
  [ vdlfl1z*] ;
  [ vdwdimmz*] ;
  [ vdwddelz* ] ;
  [ vdexf7z* ] ;
  
  MODEL SYNTHETIC:
  
  gcp BY vdori* (l1);
  gcp BY vdlfl1z (l2);
  gcp BY vdlfl2 (l3); 
  gcp BY vdlfl3 (l4);
  gcp BY vdwdimmz (l5);
  gcp BY vdwddelz (l6);
  gcp BY vdexf7z (l7);
  gcp BY vdsevens (l8);
  gcp BY vdcount (l9);
  gcp@1;
  
  recall BY vdwdimmz* (ml1);
  recall BY vdwddelz* (ml1);
  recall@1;
  
  gcp with recall@0;

  
  [ vdori$1*] (t1);
  [ vdori$2*] (t2);
  [ vdori$3*] (t3);
  [ vdori$4*] (t4);
  [ vdlfl2$1*] (t5);
  [ vdlfl2$2*] (t6);
  [ vdlfl3$1*] (t7);
  [ vdlfl3$2*] (t8);
  [ vdsevens$1*] (t9);
  [ vdsevens$2*] (t10);
  [ vdsevens$3*] (t11);
  [ vdsevens$4*] (t12);
  [ vdsevens$5*] (t13);
  [ vdcount$1*] (t14);
  
  [ vdlfl1z*] ;
  [ vdwdimmz*] ;
  [ vdwddelz* ] ;
  [ vdexf7z* ] ;
  
  [gcp@0];
  [recall@0];
  vdori@1;
  vdlfl2@1;
  vdlfl3@1;
  vdsevens@1;
  vdcount@1;

  

  ",
  usevariables = colnames(cog_merged),
  VARIABLE = "CATEGORICAL = vdori vdlfl2 vdlfl3
  vdsevens vdcount;
  GROUP = group (1 = Observed 2 = Synthetic);",
  OUTPUT = "TECH4 STDYX SVALUES;",
  ANALYSIS = "ESTIMATOR = WLSMV;
  PARAMETERIZATION = THETA;",
  rdata = cog_merged)

metric <- mplusModeler(metric_model, 
                           modelout = "metric.inp", run = TRUE)
