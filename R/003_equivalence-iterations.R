library(furrr)
library(progressr)
library(dplyr)
library(purrr)
library(MplusAutomation)
library(synthpop)

# choose a plan:
# - Windows/macOS: multisession (separate R sessions)
# - Linux/macOS (no RStudio parallel conflicts): multicore is usually fastest
future::plan(multisession, workers = max(1, future::availableCores() - 1))

handlers(global = TRUE)

run_one_iter <- function(i,
                         data,
                         base_seed = 904,
                         minnumlevels = 5) {
  # 1) synthesize (unique seed per iter)
  syn_obj <- try(
    synthpop::syn(data, seed = base_seed + i, minnumlevels = minnumlevels),
    silent = TRUE
  )
  if (inherits(syn_obj, "try-error") || is.null(syn_obj)) {
    return(tibble(iter = i,
                  paramHeader = character(),
                  param = character(),
                  low2.5 = double(),
                  est = double(),
                  up2.5 = double()))
  }
  syn_df <- syn_obj$syn
  
  # 2) build Mplus model for this dataset
  singledomainmodel_syn <- mplusObject(
    MODEL = "
      gcp BY vdori* vdlfl1z vdlfl2 vdlfl3 vdwdimmz vdwddelz vdexf7z vdsevens vdcount;
      gcp@1;

      recall BY vdwdimmz* (1);
      recall BY vdwddelz* (1);
      recall@1;

      gcp with recall@0;
    ",
    usevariables = colnames(syn_df),
    VARIABLE = "CATEGORICAL = vdori vdlfl2 vdlfl3 vdsevens vdcount;",
    OUTPUT = "TECH4 STDYX CINTERVAL MODINDICES(ALL);",
    ANALYSIS = "ESTIMATOR = WLSMV; PARAMETERIZATION = THETA;",
    rdata = syn_df
  )
  
  # 3) unique working dir + model file per iter (prevents clashes in parallel)
  out_dir   <- file.path(getwd(), "mplus_runs", sprintf("iter_%04d", i))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  inp_file  <- file.path(out_dir, sprintf("gcpsyn_%04d.inp", i))
  
  # 4) run Mplus & extract CI table (robust to failures)
  mres <- try(
    mplusModeler(
      singledomainmodel_syn,
      modelout = inp_file,
      run = TRUE
      # , writeData = "ifmissing"  # uncomment if you want to avoid re-writing data
    ),
    silent = TRUE
  )
  
  if (inherits(mres, "try-error") || is.null(mres) ||
      is.null(mres$results$parameters$ci.stdyx.standardized)) {
    return(tibble(iter = i,
                  paramHeader = character(),
                  param = character(),
                  low2.5 = double(),
                  est = double(),
                  up2.5 = double()))
  }
  
  mres$results$parameters$ci.stdyx.standardized %>%
    select(paramHeader, param, low2.5, est, up2.5) %>%
    mutate(iter = i, .before = 1)
}

n_iter <- 1001

with_progress({
  p <- progressor(steps = n_iter)
  syn_values_all <- future_map_dfr(
    1:n_iter,
    ~{
      res <- run_one_iter(
        i = .x,
        data = demo_and_cog,   # <<< your observed dataset
        base_seed = 904,
        minnumlevels = 5
      )
      p(sprintf("finished %d/%d", .x, n_iter))
      res
    },
    .options = furrr_options(seed = TRUE)  # ensure reproducible RNG in parallel
  )
})

# inspect
dplyr::glimpse(syn_values_all)

combined_values_all <- equivalence_by_iter(obs_values, syn_values_all)