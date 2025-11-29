#!/usr/bin/env Rscript
# Fit hierarchical models predicting choice from listenBPM and responseBPM interaction
# Model: Decision ~ listenBPM * responseBPM (interoceptive trials only)
# Includes maximal random slopes including the interaction term

# Ensure required packages are installed
if (!require("pacman")) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(dplyr, readr, glmmTMB, ggplot2, sjPlot, here, tidyr, emmeans)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "revision1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
DOCS_DIR <- file.path(OUTPUT_DIR, "docs")

# Create directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", OUTPUT_DIR)

# Load data
INPUT_CSV <- file.path(BASE_DIR, "data", "HRD_trial_level_data.csv")
if (!file.exists(INPUT_CSV)) stop("Input CSV not found: ", INPUT_CSV)

message("Reading data from: ", INPUT_CSV)
df <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)

# Load HRV data
HRV_CSV <- file.path(BASE_DIR, "data", "alldrugs_full_HRV.csv")
if (!file.exists(HRV_CSV)) stop("HRV CSV not found: ", HRV_CSV)

message("Reading HRV data from: ", HRV_CSV)
df_hrv <- readr::read_csv(HRV_CSV, show_col_types = FALSE)

# Prepare HRV data for merging
# Note: SubNo in HRV data (e.g., "0007") seems to correspond to the last 2 digits of the subject ID
# in the trial data (e.g., "sub_3007" or "sub_4008").
# There is some ambiguity for IDs that exist in both cohorts (e.g., 15, 24),
# but we will map based on the last 2 digits.
df_hrv_clean <- df_hrv %>%
  mutate(
    # Extract last 2 digits for matching
    numeric_id = stringr::str_sub(SubNo, -2),
    # Map Drug names to match trial data
    drugs = case_when(
      Drug == "plac" ~ "PLACEBO",
      Drug == "biso" ~ "BISO",
      Drug == "prop" ~ "PROP",
      TRUE ~ toupper(Drug)
    )
  ) %>%
  select(numeric_id, drugs, meanHR)

# 0. Check required columns
required_cols <- c("Modality", "Decision", "listenBPM", "responseBPM",
                   "subject", "drugs")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in CSV: ",
       paste(missing_cols, collapse = ", "))
}

# 1. Filter to interoceptive trials only
# Check if Modality column exists and has Intero
if ("Modality" %in% names(df)) {
  df_intero <- df %>%
    filter(Modality == "Intero")
} else {
  warning("Modality column not found, assuming all data is relevant or checking other columns.")
  df_intero <- df
}

# 2. Strict Response check
# Decision is the response column in this dataset
df_intero <- df_intero %>%
  filter(Decision %in% c("More", "Less")) %>%
  filter(!is.na(listenBPM), !is.na(responseBPM))

message("Interoceptive trials with complete data: ", nrow(df_intero))

# Convert Decision to binary (More = 1, Less = 0)
df_intero <- df_intero %>%
  mutate(Response_binary = ifelse(Decision == "More", 1, 0))

# Define groups based on drugs
# We are now combining all data, so we don't loop over drugs.
# But we need to ensure 'drugs' is a factor with PLACEBO as reference.
df_intero <- df_intero %>%
  filter(drugs %in% c("PLACEBO", "BISO", "PROP")) %>%
  mutate(drugs = factor(drugs, levels = c("PLACEBO", "BISO", "PROP")))

# Merge with HRV data
message("Merging with HRV data...")

# Create numeric_id in trial data
df_intero <- df_intero %>%
  mutate(numeric_id = stringr::str_sub(subject, -2)) %>%
  left_join(df_hrv_clean, by = c("numeric_id", "drugs"))

# Check for missing meanHR
n_missing_hr <- sum(is.na(df_intero$meanHR))
if (n_missing_hr > 0) {
  warning("Missing meanHR for ", n_missing_hr, " trials. Dropping these trials.")
  df_intero <- df_intero %>% filter(!is.na(meanHR))
}

# Re-factor drugs to ensure PLACEBO is reference (join might have reset it)
df_intero <- df_intero %>%
  mutate(drugs = factor(drugs, levels = c("PLACEBO", "BISO", "PROP")))

# Helper function for safe scaling
safe_scale <- function(x) {
  if(sd(x, na.rm=TRUE) == 0) return(rep(0, length(x)))
  return(as.vector(scale(x)))
}

message("Fitting model on combined data with Drug interaction and meanHR control...")

# Scale predictors within subjects safely
# Note: meanHR is constant within subject-drug block, but varies between subjects/drugs.
# We should probably scale it across the whole dataset or within subject?
# Since it's a control for the drug effect (which varies within subject), scaling it globally makes sense
# to interpret it as "effect of 1 SD increase in HR".
# However, listenBPM and responseBPM are trial-level.

df_intero <- df_intero %>%
  mutate(
    meanHR_scaled = safe_scale(meanHR)
  ) %>%
  group_by(subject) %>%
  mutate(
    responseBPM_scaled = safe_scale(responseBPM),
    listenBPM_scaled = safe_scale(listenBPM)
  ) %>%
  ungroup()

# Define a helper to check if a model is "good enough"
is_converged <- function(m) {
  if (inherits(m, "try-error")) return(FALSE)
  if (m$fit$convergence != 0) return(FALSE)
  if (!m$sdr$pdHess) return(FALSE)
  return(TRUE)
}

# Fit model
# Fixed: drugs * listenBPM + responseBPM + meanHR
# Random: (1 + listenBPM + responseBPM | subject)

message("Fitting maximal model...")

mod <- try(glmmTMB(
  Response_binary ~ drugs * listenBPM_scaled + responseBPM_scaled + meanHR_scaled +
    (1 + listenBPM_scaled + responseBPM_scaled | subject),
  data = df_intero,
  family = binomial(link = "logit")
), silent = TRUE)

# 2. If failed or non-converged, try Uncorrelated Random Effects
if (!is_converged(mod)) {
  message("Maximal model failed or did not converge. Trying uncorrelated random effects (||)...")
  mod <- try(glmmTMB(
    Response_binary ~ drugs * listenBPM_scaled + responseBPM_scaled + meanHR_scaled +
      (1 + listenBPM_scaled + responseBPM_scaled || subject),
    data = df_intero,
    family = binomial(link = "logit")
  ), silent = TRUE)
}

# 3. If that still fails, try random intercepts only
if (!is_converged(mod)) {
  message("Random slopes model failed. Dropping random slopes (intercept only)...")
  mod <- try(glmmTMB(
    Response_binary ~ drugs * listenBPM_scaled + responseBPM_scaled + meanHR_scaled +
      (1 | subject),
    data = df_intero,
    family = binomial(link = "logit")
  ), silent = TRUE)
}

# Final check
if (!is_converged(mod)) {
  message("WARNING: Model failed to converge fully.")
} else {
  message("Model converged successfully.")
}

# Extract summary
mod_sum <- summary(mod)

message("\nFixed effects:")
print(mod_sum$coefficients$cond)

message("\nAIC: ", AIC(mod))

# Save model object
model_file <- file.path(DATA_DIR, "model_combined_drug_interaction_with_HR.rds")
saveRDS(mod, model_file)
message("Saved model to: ", model_file)

# Save full model summary to text file
summary_text_file <- file.path(DOCS_DIR, "model_summary_combined_drug_interaction_with_HR.txt")
sink(summary_text_file)
cat("========================================\n")
cat("Model: Combined Data Drug * ListenBPM Interaction + meanHR Control\n")
cat("Formula: Response ~ drugs * listenBPM + responseBPM + meanHR + (1 + listenBPM + responseBPM | subject)\n")
cat("Subjects: ", length(unique(df_intero$subject)), "\n")
cat("Trials: ", nrow(df_intero), "\n")
cat("========================================\n\n")
print(mod_sum)
cat("\n\nAIC: ", AIC(mod), "\n")
sink()
message("Saved summary to: ", summary_text_file)

# Create summary table
coef_df <- as.data.frame(mod_sum$coefficients$cond)
coef_df$term <- rownames(coef_df)
coef_df$model <- "Combined_Drug_Interaction_HR_Control"
coef_df$AIC <- AIC(mod)

# Write summary table
summary_file <- file.path(DATA_DIR, "model_coefficients_combined_with_HR.csv")
readr::write_csv(coef_df, summary_file)
message("\nWrote summary table to: ", summary_file)

# --- Post-hoc Tests (Emmeans) ---
message("\nCalculating simple slopes for listenBPM by Drug...")

# Calculate simple slopes of listenBPM_scaled for each drug level
# We want to see if the slope of listenBPM is significant within each drug condition
emm_trends <- emmeans::emtrends(mod, specs = "drugs", var = "listenBPM_scaled")
emm_summary <- summary(emm_trends)

print(emm_summary)

# Save emmeans results
emm_file <- file.path(DATA_DIR, "emmeans_listenBPM_slopes.csv")
readr::write_csv(as.data.frame(emm_summary), emm_file)
message("Saved emmeans results to: ", emm_file)

# Pairwise comparisons of slopes
emm_pairs <- pairs(emm_trends)
emm_pairs_summary <- summary(emm_pairs)
print(emm_pairs_summary)

emm_pairs_file <- file.path(DATA_DIR, "emmeans_listenBPM_slopes_pairwise.csv")
readr::write_csv(as.data.frame(emm_pairs_summary), emm_pairs_file)
message("Saved pairwise slope comparisons to: ", emm_pairs_file)

message("\nAnalysis complete!")
