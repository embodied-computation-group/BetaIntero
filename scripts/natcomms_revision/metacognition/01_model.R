#!/usr/bin/env Rscript
# Fit ordered beta regression models of metacognitive confidence (RRST + HRD)
# RRST: Conf ~ Stimulus + Drug * Accuracy + (1 + Stimulus + Accuracy | subject)
# HRD:  Conf ~ Drug * ResponseCorrect + BPM_scaled + (1 + ResponseCorrect + BPM_scaled | subject)
# Both use family = ordbeta()

# Ensure required packages are installed
if (!require("pacman")) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(dplyr, readr, tidyr, glmmTMB, emmeans, here)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision", "metacognition")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
DOCS_DIR <- file.path(OUTPUT_DIR, "docs")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

# Create directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# Generate Timestamp
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Output directory: ", OUTPUT_DIR)
message("Run Timestamp: ", TIMESTAMP)
message("glmmTMB version: ", packageVersion("glmmTMB"))

# --- Helper: clean_and_report ---
clean_and_report <- function(data, dataset_name) {
  total_rows_before <- nrow(data)
  cleaned_data <- data %>% tidyr::drop_na()
  total_rows_after <- nrow(cleaned_data)
  dropped_pct <- ((total_rows_before - total_rows_after) / total_rows_before) * 100
  message(sprintf("Dataset: %s | Before: %d | After: %d | Dropped: %.2f%%",
                  dataset_name, total_rows_before, total_rows_after, dropped_pct))
  return(cleaned_data)
}

# ============================================================
# RRST (Respiratory) Confidence Model
# ============================================================
message("\n--- RRST Confidence Model ---")

RRST_INPUT <- file.path(BASE_DIR, "data", "cleaned", "RRST.csv")
if (!file.exists(RRST_INPUT)) stop("RRST input not found: ", RRST_INPUT)

message("Reading RRST data from: ", RRST_INPUT)
RRST_trial_data <- readr::read_csv(RRST_INPUT, show_col_types = FALSE) %>%
  filter(ConfResp >= 0 & ConfResp <= 100) %>%
  mutate(
    Conf = ConfResp / 100,
    Stimulus = Stim,
    Drug = factor(drugs,
                  levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    Condition = drugs,
    Accuracy = factor(Resp, levels = c(1, 0), labels = c("Correct", "Incorrect"))
  )

# Clean
RRST_trial_data <- clean_and_report(RRST_trial_data, "RRST_trial_data")

message("RRST subjects: ", length(unique(RRST_trial_data$subject)))
message("RRST trials: ", nrow(RRST_trial_data))

# Fit model
message("Fitting RRST ordered beta regression model...")
rrst_conf_model <- glmmTMB(
  Conf ~ Stimulus + Drug * Accuracy + (1 + Stimulus + Accuracy | subject),
  data = RRST_trial_data,
  family = ordbeta(),
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)

rrst_sum <- summary(rrst_conf_model)

message("\nRRST Fixed effects:")
print(rrst_sum$coefficients$cond)
message("\nRRST AIC: ", AIC(rrst_conf_model))

# Save model
rrst_model_filename <- sprintf("model_rrst_confidence_%s.rds", TIMESTAMP)
rrst_model_file <- file.path(DATA_DIR, rrst_model_filename)
saveRDS(rrst_conf_model, rrst_model_file)
message("Saved RRST model to: ", rrst_model_file)

# Save summary text
rrst_summary_filename <- sprintf("model_summary_rrst_confidence_%s.txt", TIMESTAMP)
rrst_summary_file <- file.path(DOCS_DIR, rrst_summary_filename)
sink(rrst_summary_file)
cat("========================================\n")
cat("Model: RRST Metacognitive Confidence (Ordered Beta Regression)\n")
cat("Run Timestamp: ", TIMESTAMP, "\n")
cat("Formula: Conf ~ Stimulus + Drug * Accuracy + (1 + Stimulus + Accuracy | subject)\n")
cat("Family: ordbeta()\n")
cat("Subjects: ", length(unique(RRST_trial_data$subject)), "\n")
cat("Trials: ", nrow(RRST_trial_data), "\n")
cat("========================================\n\n")
print(rrst_sum)
cat("\n\nAIC: ", AIC(rrst_conf_model), "\n")
sink()
message("Saved RRST summary to: ", rrst_summary_file)

# ============================================================
# HRD (Cardiac) Confidence Model
# ============================================================
message("\n--- HRD Confidence Model ---")

HRD_INPUT <- file.path(BASE_DIR, "data", "cleaned", "HRD.csv")
if (!file.exists(HRD_INPUT)) stop("HRD input not found: ", HRD_INPUT)

message("Reading HRD data from: ", HRD_INPUT)
HRD_trl_data <- read.csv(HRD_INPUT) %>%
  mutate(
    Conf = Confidence / 100,
    drugs = as.factor(drugs),
    Drug = factor(drugs,
                  levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    ResponseCorrect = factor(ResponseCorrect,
                             levels = c("True", "False"),
                             labels = c("Correct", "Incorrect")),
    BPM_scaled = scale(listenBPM),
    Condition = drugs
  ) %>%
  filter(subject != "sub_4049")

# Clean
HRD_trl_data <- clean_and_report(HRD_trl_data, "HRD_trl_data")

message("HRD subjects: ", length(unique(HRD_trl_data$subject)))
message("HRD trials: ", nrow(HRD_trl_data))

# Fit model
message("Fitting HRD ordered beta regression model...")
hrd_conf_model <- glmmTMB(
  Conf ~ Drug * ResponseCorrect + BPM_scaled + (1 + ResponseCorrect + BPM_scaled | subject),
  data = HRD_trl_data,
  family = ordbeta(),
  start = list(psi = c(0, 1)),
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)

hrd_sum <- summary(hrd_conf_model)

message("\nHRD Fixed effects:")
print(hrd_sum$coefficients$cond)
message("\nHRD AIC: ", AIC(hrd_conf_model))

# Save model
hrd_model_filename <- sprintf("model_hrd_confidence_%s.rds", TIMESTAMP)
hrd_model_file <- file.path(DATA_DIR, hrd_model_filename)
saveRDS(hrd_conf_model, hrd_model_file)
message("Saved HRD model to: ", hrd_model_file)

# Save summary text
hrd_summary_filename <- sprintf("model_summary_hrd_confidence_%s.txt", TIMESTAMP)
hrd_summary_file <- file.path(DOCS_DIR, hrd_summary_filename)
sink(hrd_summary_file)
cat("========================================\n")
cat("Model: HRD Metacognitive Confidence (Ordered Beta Regression)\n")
cat("Run Timestamp: ", TIMESTAMP, "\n")
cat("Formula: Conf ~ Drug * ResponseCorrect + BPM_scaled + (1 + ResponseCorrect + BPM_scaled | subject)\n")
cat("Family: ordbeta()\n")
cat("Subjects: ", length(unique(HRD_trl_data$subject)), "\n")
cat("Trials: ", nrow(HRD_trl_data), "\n")
cat("Outlier removed: sub_4049\n")
cat("========================================\n\n")
print(hrd_sum)
cat("\n\nAIC: ", AIC(hrd_conf_model), "\n")
sink()
message("Saved HRD summary to: ", hrd_summary_file)

# ============================================================
# Post-hoc Emmeans
# ============================================================
message("\n--- Post-hoc Emmeans ---")

# HRD: Drug * ResponseCorrect (from Analysis.Rmd)
message("HRD: emmeans for Drug * ResponseCorrect...")
hrd_emm <- emmeans(hrd_conf_model, ~ Drug * ResponseCorrect)
hrd_posthoc <- pairs(hrd_emm, by = "ResponseCorrect", adjust = "holm")

message("HRD pairwise comparisons (by ResponseCorrect, Holm-adjusted):")
print(summary(hrd_posthoc))

# Save HRD emmeans
hrd_emm_filename <- sprintf("emmeans_hrd_drug_accuracy_%s.csv", TIMESTAMP)
readr::write_csv(as.data.frame(hrd_emm), file.path(DATA_DIR, hrd_emm_filename))
message("Saved HRD emmeans to: ", hrd_emm_filename)

hrd_pairs_filename <- sprintf("emmeans_hrd_pairwise_by_accuracy_%s.csv", TIMESTAMP)
readr::write_csv(as.data.frame(summary(hrd_posthoc)), file.path(DATA_DIR, hrd_pairs_filename))
message("Saved HRD pairwise to: ", hrd_pairs_filename)

# RRST: Drug | Accuracy (from plot_3 in plots.R)
message("RRST: emmeans for Drug | Accuracy...")
rrst_emm <- emmeans(rrst_conf_model, ~ Drug | Accuracy)
rrst_posthoc <- contrast(rrst_emm, method = "pairwise")

message("RRST pairwise comparisons (by Accuracy):")
print(summary(rrst_posthoc))

# Save RRST emmeans
rrst_emm_filename <- sprintf("emmeans_rrst_drug_accuracy_%s.csv", TIMESTAMP)
readr::write_csv(as.data.frame(rrst_emm), file.path(DATA_DIR, rrst_emm_filename))
message("Saved RRST emmeans to: ", rrst_emm_filename)

rrst_pairs_filename <- sprintf("emmeans_rrst_pairwise_by_accuracy_%s.csv", TIMESTAMP)
readr::write_csv(as.data.frame(summary(rrst_posthoc)), file.path(DATA_DIR, rrst_pairs_filename))
message("Saved RRST pairwise to: ", rrst_pairs_filename)

# ============================================================
# Coefficient CSVs
# ============================================================
message("\n--- Saving coefficient tables ---")

# RRST coefficients
rrst_coef_df <- as.data.frame(rrst_sum$coefficients$cond)
rrst_coef_df$term <- rownames(rrst_coef_df)
rrst_coef_df$model <- "RRST_Confidence"
rrst_coef_df$AIC <- AIC(rrst_conf_model)
rrst_coef_df$timestamp <- TIMESTAMP

rrst_coef_filename <- sprintf("model_coefficients_rrst_confidence_%s.csv", TIMESTAMP)
readr::write_csv(rrst_coef_df, file.path(DATA_DIR, rrst_coef_filename))
message("Saved RRST coefficients to: ", rrst_coef_filename)

# HRD coefficients
hrd_coef_df <- as.data.frame(hrd_sum$coefficients$cond)
hrd_coef_df$term <- rownames(hrd_coef_df)
hrd_coef_df$model <- "HRD_Confidence"
hrd_coef_df$AIC <- AIC(hrd_conf_model)
hrd_coef_df$timestamp <- TIMESTAMP

hrd_coef_filename <- sprintf("model_coefficients_hrd_confidence_%s.csv", TIMESTAMP)
readr::write_csv(hrd_coef_df, file.path(DATA_DIR, hrd_coef_filename))
message("Saved HRD coefficients to: ", hrd_coef_filename)

message("\nAnalysis complete!")
