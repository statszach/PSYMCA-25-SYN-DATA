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

demo_and_cog <- HCAP2016 %>% 
  dplyr::filter(age65up == 1) %>% 
  dplyr::select(vdori, vdlfl1z,
                vdlfl2, vdlfl3,
                vdwdimmz, vdwddelz,
                vdexf7z, vdsevens,
                vdcount,
                rage, female, black, hisp, SCHLYRS)

# Set wd

fs::dir_create("300_random")
setwd(here::here("300_random"))

# Generate synthetic datasets

# Generate synthetic datasets with random column order each time
gen_synth_datasets <- function(n_iter = 1001, data, base_seed = 904, minnumlevels = 5) {
  synth_list <- vector("list", n_iter)
  vars <- colnames(data)   # grab variable names once
  
  for (i in seq_len(n_iter)) {
    # shuffle variable order with a fresh seed each loop
    set.seed(base_seed + i)    # ensures reproducible but different each iteration
    vars_random <- sample(vars)
    data_shuffled <- data[, vars_random]
    
    # generate synthetic dataset
    synth_list[[i]] <- synthpop::syn(
      data_shuffled,
      seed = base_seed + i,
      minnumlevels = minnumlevels
    )
  }
  
  synth_list
}

# Example call
synth_datasets <- gen_synth_datasets(
  n_iter = 1001,
  data = demo_and_cog,
  base_seed = 904,
  minnumlevels = 5
)


### Observed data

usevar <- c("vdori", "vdlfl1z", "vdlfl2", "vdlfl3", "vdwdimmz", "vdwddelz", "vdexf7z", "vdsevens", "vdcount")

singledomainmodel_obs <- mplusObject(
  MODEL = "
  gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
  gcp@1;
  
  recall BY vdwdimmz* (1);
  recall BY vdwddelz* (1);
  recall@1;
  
  gcp with recall@0;
  
  ",
  usevariables = usevar,
  VARIABLE = "CATEGORICAL = vdori vdlfl2 vdlfl3
  vdsevens vdcount;",
  OUTPUT = "TECH4 STDYX CINTERVAL SVALUES;",
  ANALYSIS = "ESTIMATOR = WLSMV;
  PARAMETERIZATION = THETA;",
  rdata = demo_and_cog)

gcp_obs <- mplusModeler(singledomainmodel_obs, 
                        modelout = "gcpobs.inp", run = TRUE)

obs_values <- gcp_obs$results$parameters$stdyx.standardized %>% 
  dplyr::select(paramHeader, param, est) %>% 
  dplyr::rename(obs_est = est)

mplus_results <- vector("list", length(synth_datasets))

for (i in seq_along(synth_datasets)) {
  
  synmplus_i <- synth_datasets[[i]]$syn  # Extract the synthetic data frame
  
  singledomainmodel_syn <- mplusObject(
    MODEL = "
      gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
      gcp@1;
      
      recall BY vdwdimmz* (1);
      recall BY vdwddelz* (1);
      recall@1;
      
      gcp with recall@0;
      
    ",
    usevariables = usevar,
    VARIABLE = "CATEGORICAL = vdori vdlfl2 vdlfl3 vdsevens vdcount;",
    OUTPUT = "TECH4 STDYX CINTERVAL MODINDICES(ALL);",
    ANALYSIS = "ESTIMATOR = WLSMV; PARAMETERIZATION = THETA;",
    rdata = synmplus_i
  )
  
  # Run Mplus and store the result
  mplus_results[[i]] <- mplusModeler(
    singledomainmodel_syn,
    modelout = paste0("gcpsyn_", i, ".inp"),  # unique file for each run
    run = TRUE
  )
}

extract_ci_one <- function(mres, iter) {
  out <- try({
    mres$results$parameters$ci.stdyx.standardized %>%
      select(paramHeader, param, low2.5, est, up2.5) %>%
      mutate(iter = iter, .before = 1)
  }, silent = TRUE)
  
  if (inherits(out, "try-error") || is.null(out)) {
    # return a zero-row tibble with expected columns (so bind_rows works)
    tibble(iter = iter,
           paramHeader = character(),
           param = character(),
           low2.5 = double(),
           est = double(),
           up2.5 = double())
  } else {
    out
  }
}

# Apply to every Mplus run and bind
syn_values_all <- map2_dfr(mplus_results, seq_along(mplus_results), extract_ci_one)

equivalence_by_iter <- function(obs_values, syn_values_all,
                                cal_strict = c(0.05, 0.10),
                                cal_accept = c(0.10, 0.20)) {
  iters <- sort(unique(syn_values_all$iter))
  
  # replicate obs rows for every iteration
  base <- tidyr::crossing(iter = iters, obs_values)
  
  # join to the synthetic CI table (keeps iter-specific rows)
  combined <- base %>%
    left_join(
      syn_values_all %>%
        dplyr::select(iter, paramHeader, param, low2.5, est, up2.5),
      by = c("iter", "paramHeader", "param")
    ) %>%
    # your equivalence rules
    mutate(
      strict_equivalence = dplyr::case_when(
        paramHeader %in% c("GCP.BY", "RECALL.BY", "Residual.Variances") &
          (low2.5 >= obs_est - cal_strict[1] & up2.5 <= obs_est + cal_strict[1]) ~ "PASS",
        paramHeader %in% c("Thresholds", "Intercepts") &
          (low2.5 >= obs_est - cal_strict[2] & up2.5 <= obs_est + cal_strict[2]) ~ "PASS",
        TRUE ~ "FAIL"
      ),
      acceptable_equivalence = dplyr::case_when(
        paramHeader %in% c("GCP.BY", "RECALL.BY", "Residual.Variances") &
          (low2.5 >= obs_est - cal_accept[1] & up2.5 <= obs_est + cal_accept[1]) ~ "PASS",
        paramHeader %in% c("Thresholds", "Intercepts") &
          (low2.5 >= obs_est - cal_accept[2] & up2.5 <= obs_est + cal_accept[2]) ~ "PASS",
        TRUE ~ "FAIL"
      )
    )
  
  combined
}

# Run it
combined_values_all <- equivalence_by_iter(obs_values, syn_values_all)

equiv_summary <- combined_values_all %>%
  group_by(paramHeader, param) %>%
  summarise(
    strict_pass_rate = mean(strict_equivalence == "PASS", na.rm = TRUE),
    acceptable_pass_rate = mean(acceptable_equivalence == "PASS", na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )