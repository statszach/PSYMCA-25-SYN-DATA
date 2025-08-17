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
  gtsummary::add_difference(list(all_continuous() ~ "cohens_d")) %>% 
  gtsummary::as_gt() %>% 
  gt::gtsave("t1.rtf")

#### Sarah will show how to do equivalence testing
#### to get a better idea of how to tell if
#### there are differences in synthetic vs real




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
  OUTPUT = "TECH4 STDYX CINTERVAL SVALUES;",
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

## Build table to check equivalence

obs_values <- gcp_obs$results$parameters$stdyx.standardized %>% 
  dplyr::select(paramHeader, param, est) %>% 
  dplyr::rename(obs_est = est)

syn_values <- gcp_syn$results$parameters$ci.stdyx.standardized %>% 
  dplyr::select(paramHeader, param, low2.5, est, up2.5)

combined_values <- obs_values %>% 
  left_join(syn_values, by = c("paramHeader", "param")) %>% 
  dplyr::mutate(strict_equivalence = dplyr::case_when(
    paramHeader %in% c("GCP.BY", "RECALL.BY", "Residual.Variances") &
      (low2.5 >= obs_est - 0.05 & up2.5 <= obs_est + 0.05) ~ "PASS",
    
    paramHeader %in% c("Thresholds", "Intercepts") &
      (low2.5 >= obs_est - 0.10 & up2.5 <= obs_est + 0.10) ~ "PASS",
    
    TRUE ~ "FAIL"),
      acceptable_equivalence = dplyr::case_when(
        paramHeader %in% c("GCP.BY", "RECALL.BY", "Residual.Variances") &
          (low2.5 >= obs_est - 0.1 & up2.5 <= obs_est + 0.1) ~ "PASS",
        
        paramHeader %in% c("Thresholds", "Intercepts") &
          (low2.5 >= obs_est - 0.2 & up2.5 <= obs_est + 0.2) ~ "PASS",
        
        TRUE ~ "FAIL"))

# For loadings, use .05 and .1
# For thresholds, use .10 and .20
# For residual variances, use ??



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
  rdata = cog_merged,
  SAVEDATA = "H5RESULTS = results;
              DIFFTEST = configural.dat;")

configural <- mplusModeler(configural_model, 
                        modelout = "configural.inp", run = TRUE)

fullyconstrained_model <- mplusObject(
  MODEL = "
  
  gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
  gcp@1;
  
  recall BY vdwdimmz* (1);
  recall BY vdwddelz* (1);
  recall@1;
  
  gcp with recall@0;
  
  MODEL OBSERVED:
  
  gcp BY vdori*0.64973 (2);
  gcp BY vdlfl1z*0.04753 (3);
  gcp BY vdlfl2*1.05251 (4);
  gcp BY vdlfl3*0.84640 (5);
  gcp BY vdwdimmz*0.10688 (6);
  gcp BY vdwddelz*0.12884 (7);
  gcp BY vdexf7z*0.11547 (8);
  gcp BY vdsevens*0.90846 (9);
  gcp BY vdcount*0.72949 (10);
  recall BY vdwdimmz*0.11833 (33);
  recall BY vdwddelz*0.11833 (33);

  gcp WITH recall@0;

  [ vdlfl1z*0.16760 ] (11);
  [ vdwdimmz*0.49698 ] (12);
  [ vdwddelz*0.39135 ] (13);
  [ vdexf7z*0.60927 ] (14);
  [ gcp@0 ];
  [ recall@0 ];

  [ vdori$1*-2.63384 ] (15);
  [ vdori$2*-2.25207 ] (16);
  [ vdori$3*-1.81333 ] (17);
  [ vdori$4*-0.77681 ] (18);
  [ vdlfl2$1*-3.55421 ] (19);
  [ vdlfl2$2*-1.93039 ] (20);
  [ vdlfl3$1*-2.39380 ] (21);
  [ vdlfl3$2*-0.50783 ] (22);
  [ vdsevens$1*-1.65874 ] (23);
  [ vdsevens$2*-1.03448 ] (24);
  [ vdsevens$3*-0.67764 ] (25);
  [ vdsevens$4*-0.22317 ] (26);
  [ vdsevens$5*0.40499 ] (27);
  [ vdcount$1*-1.72961 ] (28);

  vdori@1;
  vdlfl1z*0.00360 (29);
  vdlfl2@1;
  vdlfl3@1;
  vdwdimmz*0.00471 (30);
  vdwddelz*0.01249 (31);
  vdexf7z*0.02277 (32);
  vdsevens@1 ;
  vdcount@1;
  gcp@1;
  recall@1;
  
  MODEL SYNTHETIC:
  
  gcp BY vdori*0.64973 (2);
  gcp BY vdlfl1z*0.04753 (3);
  gcp BY vdlfl2*1.05251 (4);
  gcp BY vdlfl3*0.84640 (5);
  gcp BY vdwdimmz*0.10688 (6);
  gcp BY vdwddelz*0.12884 (7);
  gcp BY vdexf7z*0.11547 (8);
  gcp BY vdsevens*0.90846 (9);
  gcp BY vdcount*0.72949 (10);
  recall BY vdwdimmz*0.11833 (33);
  recall BY vdwddelz*0.11833 (33);

  gcp WITH recall@0;

  [ vdlfl1z*0.16760 ] (11);
  [ vdwdimmz*0.49698 ] (12);
  [ vdwddelz*0.39135 ] (13);
  [ vdexf7z*0.60927 ] (14);
  [ gcp@0 ];
  [ recall@0 ];

  [ vdori$1*-2.63384 ] (15);
  [ vdori$2*-2.25207 ] (16);
  [ vdori$3*-1.81333 ] (17);
  [ vdori$4*-0.77681 ] (18);
  [ vdlfl2$1*-3.55421 ] (19);
  [ vdlfl2$2*-1.93039 ] (20);
  [ vdlfl3$1*-2.39380 ] (21);
  [ vdlfl3$2*-0.50783 ] (22);
  [ vdsevens$1*-1.65874 ] (23);
  [ vdsevens$2*-1.03448 ] (24);
  [ vdsevens$3*-0.67764 ] (25);
  [ vdsevens$4*-0.22317 ] (26);
  [ vdsevens$5*0.40499 ] (27);
  [ vdcount$1*-1.72961 ] (28);

  vdori@1;
  vdlfl1z*0.00360 (29);
  vdlfl2@1;
  vdlfl3@1;
  vdwdimmz*0.00471 (30);
  vdwddelz*0.01249 (31);
  vdexf7z*0.02277 (32);
  vdsevens@1 ;
  vdcount@1;
  gcp@1;
  recall@1;

  

  ",
  usevariables = colnames(cog_merged),
  VARIABLE = "CATEGORICAL = vdori vdlfl2 vdlfl3
  vdsevens vdcount;
  GROUP = group (1 = Observed 2 = Synthetic);",
  OUTPUT = "TECH4 STDYX SVALUES;",
  ANALYSIS = "ESTIMATOR = WLSMV;
  PARAMETERIZATION = THETA;",
  rdata = cog_merged,
  SAVEDATA = "H5RESULTS = results;")

fullyconstrained <- mplusModeler(fullyconstrained_model, 
                           modelout = "fullyconstrained.inp", run = TRUE)

fullyconstraineddifftest_model <- mplusObject(
  MODEL = "
  
  gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
  gcp@1;
  
  recall BY vdwdimmz* (1);
  recall BY vdwddelz* (1);
  recall@1;
  
  gcp with recall@0;
  
  MODEL OBSERVED:
  
  gcp BY vdori*0.64973 (2);
  gcp BY vdlfl1z*0.04753 (3);
  gcp BY vdlfl2*1.05251 (4);
  gcp BY vdlfl3*0.84640 (5);
  gcp BY vdwdimmz*0.10688 (6);
  gcp BY vdwddelz*0.12884 (7);
  gcp BY vdexf7z*0.11547 (8);
  gcp BY vdsevens*0.90846 (9);
  gcp BY vdcount*0.72949 (10);
  recall BY vdwdimmz*0.11833 (33);
  recall BY vdwddelz*0.11833 (33);

  gcp WITH recall@0;

  [ vdlfl1z*0.16760 ] (11);
  [ vdwdimmz*0.49698 ] (12);
  [ vdwddelz*0.39135 ] (13);
  [ vdexf7z*0.60927 ] (14);
  [ gcp@0 ];
  [ recall@0 ];

  [ vdori$1*-2.63384 ] (15);
  [ vdori$2*-2.25207 ] (16);
  [ vdori$3*-1.81333 ] (17);
  [ vdori$4*-0.77681 ] (18);
  [ vdlfl2$1*-3.55421 ] (19);
  [ vdlfl2$2*-1.93039 ] (20);
  [ vdlfl3$1*-2.39380 ] (21);
  [ vdlfl3$2*-0.50783 ] (22);
  [ vdsevens$1*-1.65874 ] (23);
  [ vdsevens$2*-1.03448 ] (24);
  [ vdsevens$3*-0.67764 ] (25);
  [ vdsevens$4*-0.22317 ] (26);
  [ vdsevens$5*0.40499 ] (27);
  [ vdcount$1*-1.72961 ] (28);

  vdori@1;
  vdlfl1z*0.00360 (29);
  vdlfl2@1;
  vdlfl3@1;
  vdwdimmz*0.00471 (30);
  vdwddelz*0.01249 (31);
  vdexf7z*0.02277 (32);
  vdsevens@1 ;
  vdcount@1;
  gcp@1;
  recall@1;
  
  MODEL SYNTHETIC:
  
  gcp BY vdori*0.64973 (2);
  gcp BY vdlfl1z*0.04753 (3);
  gcp BY vdlfl2*1.05251 (4);
  gcp BY vdlfl3*0.84640 (5);
  gcp BY vdwdimmz*0.10688 (6);
  gcp BY vdwddelz*0.12884 (7);
  gcp BY vdexf7z*0.11547 (8);
  gcp BY vdsevens*0.90846 (9);
  gcp BY vdcount*0.72949 (10);
  recall BY vdwdimmz*0.11833 (33);
  recall BY vdwddelz*0.11833 (33);

  gcp WITH recall@0;

  [ vdlfl1z*0.16760 ] (11);
  [ vdwdimmz*0.49698 ] (12);
  [ vdwddelz*0.39135 ] (13);
  [ vdexf7z*0.60927 ] (14);
  [ gcp@0 ];
  [ recall@0 ];

  [ vdori$1*-2.63384 ] (15);
  [ vdori$2*-2.25207 ] (16);
  [ vdori$3*-1.81333 ] (17);
  [ vdori$4*-0.77681 ] (18);
  [ vdlfl2$1*-3.55421 ] (19);
  [ vdlfl2$2*-1.93039 ] (20);
  [ vdlfl3$1*-2.39380 ] (21);
  [ vdlfl3$2*-0.50783 ] (22);
  [ vdsevens$1*-1.65874 ] (23);
  [ vdsevens$2*-1.03448 ] (24);
  [ vdsevens$3*-0.67764 ] (25);
  [ vdsevens$4*-0.22317 ] (26);
  [ vdsevens$5*0.40499 ] (27);
  [ vdcount$1*-1.72961 ] (28);

  vdori@1;
  vdlfl1z*0.00360 (29);
  vdlfl2@1;
  vdlfl3@1;
  vdwdimmz*0.00471 (30);
  vdwddelz*0.01249 (31);
  vdexf7z*0.02277 (32);
  vdsevens@1 ;
  vdcount@1;
  gcp@1;
  recall@1;

  

  ",
  usevariables = colnames(cog_merged),
  VARIABLE = "CATEGORICAL = vdori vdlfl2 vdlfl3
  vdsevens vdcount;
  GROUP = group (1 = Observed 2 = Synthetic);",
  OUTPUT = "TECH4 STDYX SVALUES;",
  ANALYSIS = "ESTIMATOR = WLSMV;
  PARAMETERIZATION = THETA;
  DIFFTEST is configural.dat ; ",
  rdata = cog_merged,
  SAVEDATA = "H5RESULTS = results;")

fullyconstraineddifftest <- mplusModeler(fullyconstraineddifftest_model, 
                                 modelout = "fullyconstraineddifftest.inp", run = TRUE)

#### Test of covariances

covarfree_model <- mplusObject(
  MODEL = "
  

    MODEL OBSERVED:

    vdori vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
    
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
    
    vdori with vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
    vdlfl1z with vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
    vdlfl2 with vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
    vdlfl3 with vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
    vdwdimmz with vdwddelz vdexf7z vdsevens vdcount;
    vdwddelz with vdexf7z vdsevens vdcount;
    vdexf7z with vdsevens vdcount;
    vdsevens with vdcount;
    
    MODEL SYNTHETIC:

    vdori vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
    
    
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
    
    vdori with vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
    vdlfl1z with vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
    vdlfl2 with vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
    vdlfl3 with vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
    vdwdimmz with vdwddelz vdexf7z vdsevens vdcount;
    vdwddelz with vdexf7z vdsevens vdcount;
    vdexf7z with vdsevens vdcount;
    vdsevens with vdcount;

  

  ",
  usevariables = colnames(cog_merged),
  VARIABLE = "CATEGORICAL = vdori vdlfl2 vdlfl3
  vdsevens vdcount;
  GROUP = group (1 = Observed 2 = Synthetic);",
  OUTPUT = "TECH4 STDYX SVALUES;",
  ANALYSIS = "ESTIMATOR = WLSMV;
  PARAMETERIZATION = THETA;",
  rdata = cog_merged,
  SAVEDATA = "DIFFTEST = freecovar.dat;")

covarfree <- mplusModeler(covarfree_model, 
                           modelout = "covarfree.inp", run = TRUE)


#### Now do Rich's DIF magnitude testing

## This function works for continuous items

h5is <- here::here("results.h5")
itemis <- c("VDEXF7Z") #VDLFL1Z, "VDWDIMMZ", "VDWDDELZ", "")
group1is <- "OBSERVED"
group2is <- "SYNTHETIC"
factoris <- "GCP"
groupvaris <- "group"

params <- pull_itemParameters(h5is,group1is,group2is,itemis,factoris)
musigma <- pull_muSigma(h5is,group2is,factoris)
pooled_sd_value <- cog_merged %>%  pooled_sd(groupvaris, itemis)

results <- computeAreas(params, musigma, sd = pooled_sd_value, cov_matrix = NULL)

# vdlfl1z_d <- results$std_unsigned_area
# vdwdimmz_d <- results$std_unsigned_area
# vdwddelz_d <- results$std_unsigned_area
# vdexf7z_d <- results$std_unsigned_area

## This function works for categorical items

# Orientation

fis <- seq(-4, 4, by = 0.01) 
oria1 <- c(-2.634, -2.252, -1.813, -0.777)
orib1 <- 0.650
oria2 <- c(-2.578, -2.232, -1.800, -0.774)
orib2 <- 0.631

ori_obs <- expectedScore(fis, oria1, orib1)
ori_syn <- expectedScore(fis, oria2, orib2)

oriresults <- areaMeasures(ori_obs, ori_syn, 0, 1, fis)

vdori_h <- oriresults$SAh

# LFL2

lfl2a1 <- c(-3.554, -1.930)
lfl2b1 <- 1.053
lfl2a2 <- c(-3.221, -1.828)
lfl2b2 <- 0.902

lfl2_obs <- expectedScore(fis, lfl2a1, lfl2b1)
lfl2_syn <- expectedScore(fis, lfl2a2, lfl2b2)

lfl2results <- areaMeasures(lfl2_obs, lfl2_syn, 0, 1, fis)

vdlfl2_h <- lfl2results$SAh

# LFL3

lfl3a1 <- c(-2.394, -0.508)
lfl3b1 <- 0.846
lfl3a2 <- c(-2.312, -0.449)
lfl3b2 <- 0.791

lfl3_obs <- expectedScore(fis, lfl3a1, lfl3b1)
lfl3_syn <- expectedScore(fis, lfl3a2, lfl3b2)

lfl3results <- areaMeasures(lfl3_obs, lfl3_syn, 0, 1, fis)

vdlfl3_h <- lfl3results$SAh

# SEVENS

sevensa1 <- c(-1.659, -1.034, -0.678, -0.223, 0.405)
sevensb1 <- 0.908
sevensa2 <- c(-1.689, -1.053, -0.703, -0.226, 0.416)
sevensb2 <- 0.911

sevens_obs <- expectedScore(fis, sevensa1, sevensb1)
sevens_syn <- expectedScore(fis, sevensa2, sevensb2)

sevensresults <- areaMeasures(sevens_obs, sevens_syn, 0, 1, fis)

vdsevens_h <- sevensresults$SAh

# COUNT

counta1 <- -1.730
countb1 <- 0.729
counta2 <- -1.703
countb2 <- 0.661

count_obs <- expectedScore(fis, counta1, countb1)
count_syn <- expectedScore(fis, counta2, countb2)

countresults <- areaMeasures(count_obs, count_syn, 0, 1, fis)

vdcount_h <- countresults$SAh

# Table of results

t3 <- data.frame(
  
  item = c("VDORI",
           "VDLFL1Z",
           "VDLFL2",
           "VDLFL3",
           "VDWDIMMZ",
           "VDWDDELZ",
           "VDEXF7Z",
           "VDSEVENS",
           "VDCOUNT"),
  dorh = c(vdori_h,
           vdlfl1z_d,
           vdlfl2_h,
           vdlfl3_h,
           vdwdimmz_d,
           vdwddelz_d,
           vdexf7z_d,
           vdsevens_h,
           vdcount_h))

# Let's make a plot

plotdata <- combined_values %>% 
  dplyr::select(paramHeader, param, obs_est, low2.5, est, up2.5) %>% 
  dplyr::filter(paramHeader != "GCP.WITH") %>% 
  dplyr::filter(paramHeader != "Variances") %>% 
  dplyr::mutate(label = paste(paramHeader, param)) %>% 
  dplyr::mutate(est = dplyr::case_when(paramHeader == "GCP.BY" ~ atanh(est),
                                       paramHeader == "RECALL.BY" ~ atanh(est),
                                       TRUE ~ est),
                obs_est = dplyr::case_when(paramHeader == "GCP.BY" ~ atanh(obs_est),
                                       paramHeader == "RECALL.BY" ~ atanh(obs_est),
                                       TRUE ~ obs_est))

intthresh <- plotdata %>% 
  dplyr::filter(paramHeader == "Thresholds" | paramHeader == "Intercepts")

loadings <- plotdata %>% 
  dplyr::filter(paramHeader == "GCP.BY" | paramHeader == "RECALL.BY")

resvar <- plotdata %>% 
  dplyr::filter(paramHeader == "Residual.Variances")

ggplot(data = plotdata, aes(y = label, x = obs_est)) +
  theme_minimal() +
  geom_point(data = plotdata, aes(y = label, x = est), color = "purple", size = 2) +
  geom_point(color = "black") +
  xlab("Estimate") +
  ylab("Parameters")

ggplot(data = plotdata, aes(x = est, y = obs_est)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey50") +
  geom_point() +
  theme_minimal() +
  xlab("Synthetic Estimate") +
  ylab("Observed Estimate") +
  xlim(-4, 4) +
  ylim(-4, 4)

ggplot(data = intthresh, aes(y = label, x = obs_est)) +
  theme_minimal() +
  geom_pointrange(data = intthresh, aes(y = label,
                  x = est, xmin = low2.5, xmax = up2.5), color = "purple") +
  geom_point(color = "black") +
  xlab("Estimate") +
  ylab("Parameters")

ggplot(data = loadings, aes(y = label, x = obs_est)) +
  theme_minimal() +
  geom_pointrange(data = loadings, aes(y = label,
                                        x = est, xmin = low2.5, xmax = up2.5), color = "purple") +
  geom_point(color = "black") +
  xlab("Estimate") +
  ylab("Parameters")


ggplot(data = resvar, aes(y = label, x = obs_est)) +
  theme_minimal() +
  geom_pointrange(data = resvar, aes(y = label,
                                       x = est, xmin = low2.5, xmax = up2.5), color = "purple") +
  geom_point(color = "black") +
  xlab("Estimate") +
  ylab("Parameters")

