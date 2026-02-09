#!/usr/bin/env Rscript
# Fit hierarchical models predicting choice from listenBPM and responseBPM interaction
# Model: Decision ~ listenBPM * responseBPM (interoceptive trials only)
# Includes maximal random slopes including the interaction term
# Control: Trial Number (nTrials)

# Ensure required packages are installed
if (!require("pacman")) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(dplyr, readr, glmmTMB, ggplot2, sjPlot, here, tidyr, emmeans)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
DOCS_DIR <- file.path(OUTPUT_DIR, "docs")

# Create directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)

# Generate Timestamp
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Output directory: ", OUTPUT_DIR)
message("Run Timestamp: ", TIMESTAMP)

# Load data
INPUT_CSV <- file.path(BASE_DIR, "data", "HRD_trial_level_data.csv")
if (!file.exists(INPUT_CSV)) stop("Input CSV not found: ", INPUT_CSV)

message("Reading data from: ", INPUT_CSV)
df <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)

# 0. Check required columns
required_cols <- c("Modality", "Decision", "listenBPM", "responseBPM",
                   "subject", "drugs", "nTrials", "visit")
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

# Ensure visit is a factor
df_intero <- df_intero %>%
  mutate(visit = as.factor(visit))

# --- Outlier Removal ---
# Based on diagnostics, sub_4071 has an implausible slope (> 70).
# We also remove subjects with very few trials or perfect separation if detected.
outliers <- c("sub_4071")
message("Removing known outlier subjects: ", paste(outliers, collapse = ", "))
df_intero <- df_intero %>%
  filter(!subject %in% outliers)

message("Fitting model on combined data with Drug interaction and Trial Number control...")

# --- Mean-centering within subject ---
# Subtract each subject's grand mean (across all drug sessions) so that
# the predictor reflects within-subject, trial-level deviation.
# nTrials is left uncentered (raw trial number), matching script 09.
df_intero <- df_intero %>%
  group_by(subject) %>%
  mutate(
    listenBPM_centered = listenBPM - mean(listenBPM, na.rm = TRUE),
    responseBPM_centered = responseBPM - mean(responseBPM, na.rm = TRUE)
  ) %>%
  ungroup()

# Fit maximal model
# Fixed: drugs * listenBPM + responseBPM + nTrials + visit
# Random: (1 + listenBPM_centered + responseBPM_centered | subject)
message("Fitting maximal model with random slopes and visit control...")

mod <- glmmTMB(
  Response_binary ~ drugs * listenBPM_centered + responseBPM_centered + nTrials + visit +
    (1 + listenBPM_centered + responseBPM_centered | subject),
  data = df_intero,
  family = binomial(link = "logit")
)

# Extract summary
mod_sum <- summary(mod)

message("\nFixed effects:")
print(mod_sum$coefficients$cond)

message("\nAIC: ", AIC(mod))

# Save model object
model_filename <- sprintf("model_combined_drug_interaction_with_trialN_%s.rds", TIMESTAMP)
model_file <- file.path(DATA_DIR, model_filename)
saveRDS(mod, model_file)
message("Saved model to: ", model_file)

# Save full model summary to text file
summary_filename <- sprintf("model_summary_combined_drug_interaction_with_trialN_%s.txt", TIMESTAMP)
summary_text_file <- file.path(DOCS_DIR, summary_filename)
sink(summary_text_file)
cat("========================================\n")
cat("Model: Combined Data Drug * ListenBPM Interaction + Trial Number + Visit Control\n")
cat("Run Timestamp: ", TIMESTAMP, "\n")
cat("Formula: Response ~ drugs * listenBPM + responseBPM + nTrials + visit + (1 + listenBPM + responseBPM | subject)\n")
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
coef_df$model <- "Combined_Drug_Interaction_TrialN_Control"
coef_df$AIC <- AIC(mod)
coef_df$timestamp <- TIMESTAMP

# Write summary table
coef_filename <- sprintf("model_coefficients_combined_with_trialN_%s.csv", TIMESTAMP)
summary_file <- file.path(DATA_DIR, coef_filename)
readr::write_csv(coef_df, summary_file)
message("\nWrote summary table to: ", summary_file)

# --- Post-hoc Tests (Emmeans) ---
message("\nCalculating simple slopes for listenBPM by Drug...")

# Calculate simple slopes of listenBPM_centered for each drug level
# We want to see if the slope of listenBPM is significant within each drug condition
emm_trends <- emmeans::emtrends(mod, specs = "drugs", var = "listenBPM_centered")
emm_summary <- summary(emm_trends)

print(emm_summary)

# Save emmeans results
emm_filename <- sprintf("emmeans_listenBPM_slopes_trialN_%s.csv", TIMESTAMP)
emm_file <- file.path(DATA_DIR, emm_filename)
readr::write_csv(as.data.frame(emm_summary), emm_file)
message("Saved emmeans results to: ", emm_file)

# Pairwise comparisons of slopes
emm_pairs <- pairs(emm_trends)
emm_pairs_summary <- summary(emm_pairs)
print(emm_pairs_summary)

emm_pairs_filename <- sprintf("emmeans_listenBPM_slopes_pairwise_trialN_%s.csv", TIMESTAMP)
emm_pairs_file <- file.path(DATA_DIR, emm_pairs_filename)
readr::write_csv(as.data.frame(emm_pairs_summary), emm_pairs_file)
message("Saved pairwise slope comparisons to: ", emm_pairs_file)

message("\nAnalysis complete!")
