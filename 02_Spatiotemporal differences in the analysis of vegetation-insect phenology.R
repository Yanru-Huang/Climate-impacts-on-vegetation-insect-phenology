library(dplyr)
library(readr)
library(lme4)
library(lmerTest)
library(purrr)
library(stringr)
library(performance)
library(MuMIn)

conflict_prefer("filter", "dplyr")
conflict_prefer(lmerTest::lmer)  

# ===============================================================
#  1. Paths & Global Parameters ----------------------------------
# ===============================================================
DATA_DIR   <- ".../Coleoptera"    # Directory containing input CSV files
OUTPUT_DIR <- ".../Coleoptera2"   # Directory for output results


GRID      <- 1      
MIN_OBS   <- 30     
MIN_YEARS <- 10     
VAR_EPS   <- 1e-8   # Variance threshold for year0
CATE_VAL  <- 1      # Phenological pattern to analyse 1 or 2

# List all CSV files in the data directory
csv_files <- list.files(DATA_DIR, pattern = "\\.csv$", full.names = TRUE)


# Function to fit a Linear Mixed Model (LMM) or fallback to OLS
model_species <- function(csv_path) {
  message(">>> processing: ", csv_path)
  
  raw_df <- read_csv(csv_path, show_col_types = FALSE)
  
  ## ----------------------- Preprocessing -----------------------
  df <- raw_df %>%
    filter(cate == CATE_VAL) %>%                              # Keep only specified category
    mutate(
      across(c(Latitude, Longitude, julian_ins), as.numeric),
      year = as.integer(year)
    ) %>%
    filter(
      between(Latitude , 32, 72),                            
      between(Longitude, -15, 36),                          
      between(julian_ins, 1, 366),                          
      complete.cases(Latitude, Longitude, julian_ins, year)  
    ) %>%
    mutate(
      year0   = year - mean(year),                          
      lat_bin = floor(Latitude  / GRID) * GRID,            
      lon_bin = floor(Longitude / GRID) * GRID,            
      gridID  = paste(lat_bin, lon_bin, sep = "_")         
    )
  n_all = nrow(df)                                          
  
  # Filter grids by minimum observations and years
  grid_info <- df %>%
    group_by(gridID, lat_bin, lon_bin) %>%
    summarise(
      n_obs   = n(),                                        # Number of records in each grid
      n_years = n_distinct(year),                           # Number of distinct years
      .groups = "drop"
    ) %>%
    filter(n_obs >= MIN_OBS, n_years >= MIN_YEARS)
  
  df2 <- semi_join(df, grid_info, by = "gridID")
  
  if (nrow(df) < 400) {
    warning("No data; skipped.")
    return(invisible(NULL))
  }
  
  ## ----------------------- LMM Fitting -------------------------
  lmm_fit <- tryCatch(
    lmer(julian_ins ~ year0 + (year0 | gridID), data = df2, REML = FALSE),
    error = function(e) NULL
  )
  
  use_lmm <- !is.null(lmm_fit) && !isSingular(lmm_fit)
  
  ## ----------------------- Compile Results ---------------------
  if (use_lmm) {
    # Fixed effect estimate and p-value
    fe_tab   <- coef(summary(lmm_fit))
    fe_est   <- fe_tab["year0", "Estimate"]
    fe_pval  <- fe_tab["year0", "Pr(>|t|)"]
    
    # Random effects and p-values for slopes
    re_mat <- ranef(lmm_fit, condVar = TRUE)$gridID
    pv_arr <- attr(re_mat, "postVar")
    dimnames(pv_arr) <- list(colnames(re_mat), colnames(re_mat), NULL)
    
    rand_df <- tibble(
      gridID     = rownames(re_mat),
      rand_slope = re_mat[, "year0"],                      #random slope
      se_slope   = sqrt(map_dbl(seq_len(dim(pv_arr)[3]),    # Conditional standard error
                                ~ pv_arr["year0", "year0", .x])),
      z_rand     = rand_slope / se_slope,
      p_rand     = 2 * pnorm(-abs(z_rand))                   # Two-sided p-value
    ) %>%
      left_join(grid_info, by = "gridID") %>%
      mutate(
        model       = "LMM",
        slope_fixed = fe_est,
        p_fixed     = fe_pval
      ) %>%
      select(model, gridID, lat_bin, lon_bin,
             slope_fixed, p_fixed,
             rand_slope, p_rand, n_obs, n_years)
    
    out_df <- rand_df
    
  } else {
    # Fallback to ordinary least squares (OLS) if LMM fails
    lm_fit <- lm(julian_ins ~ year0 + Latitude + Longitude, data = df)
    co     <- summary(lm_fit)$coefficients
    
    out_df <- tibble(
      model        = "OLS",
      gridID       = NA_character_,
      lat_bin      = NA_real_,
      lon_bin      = NA_real_,
      slope_fixed  = co["year0", "Estimate"],
      p_fixed      = co["year0", "Pr(>|t|)"],
      rand_slope   = NA_real_,
      p_rand       = NA_real_,
      n_obs        = nrow(df),
      n_years      = n_distinct(df$year)
    )
  }
  
  # Add species metadata and compute actual slope
  out_df <- out_df %>%
    mutate(
      actual_slope = slope_fixed + coalesce(rand_slope, 0),
      species      = unique(raw_df$species)[1],
      order        = unique(raw_df$order)[1],
      cate         = CATE_VAL,
      n_all        = n_all
    ) %>%
    relocate(species, order, cate, .before = model)
  
  ## ----------------------- Write Output File -------------------
  sp_name   <- paste0(out_df$species[1], "_", out_df$cate[1])
  out_path  <- file.path(OUTPUT_DIR, paste0(sp_name, "_trend.csv"))
  write_csv(out_df, out_path)
  message(">>> saved: ", out_path)
  
  invisible(out_df)
}

# ===============================================================
#  3. Batch Execution & Aggregation -----------------------------
# ===============================================================
csv_files <- list.files(DATA_DIR, pattern="\\.csv$", full.names=TRUE)
results   <- purrr::compact(purrr::map(csv_files, model_species))
