#!/usr/bin/env Rscript
# Fit hierarchical model predicting choice from raw listenBPM and responseBPM
# Model: Response ~ drugs * listenBPM + responseBPM + nTrials + visit
# Random Effects: Maximal random slopes (1 + listenBPM + responseBPM | subject)
# Data: Unscaled (Raw BPM)

# Ensure required packages are installed
if (!require("pacman")) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(dplyr, readr, glmmTMB, emmeans, here)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "revision1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
DOCS_DIR <- file.path(OUTPUT_DIR, "docs")

# Create directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)

# Generate Timestamp
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Run Timestamp: ", TIMESTAMP)

# Load data
INPUT_CSV <- file.path(BASE_DIR, "data", "HRD_trial_level_data.csv")
if (!file.exists(INPUT_CSV)) stop("Input CSV not found: ", INPUT_CSV)

message("Reading data from: ", INPUT_CSV)
df <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)

# 1. Filter to interoceptive trials only
if ("Modality" %in% names(df)) {
  df_intero <- df %>% filter(Modality == "Intero")
} else {
  df_intero <- df
}

# 2. Strict Response check & Data Cleaning
df_intero <- df_intero %>%
  filter(Decision %in% c("More", "Less")) %>%
  filter(!is.na(listenBPM), !is.na(responseBPM)) %>%
  mutate(Response_binary = ifelse(Decision == "More", 1, 0)) %>%
  filter(drugs %in% c("PLACEBO", "BISO", "PROP")) %>%
  mutate(drugs = factor(drugs, levels = c("PLACEBO", "BISO", "PROP"))) %>%
  mutate(visit = as.factor(visit))

message("Interoceptive trials with complete data: ", nrow(df_intero))

# --- Outlier Removal ---
outliers <- c("sub_4071")
message("Removing known outlier subjects: ", paste(outliers, collapse = ", "))
df_intero <- df_intero %>%
  filter(!subject %in% outliers)

# --- Centering ---

# Group-mean centering: center within each subject and drug condition
message("Applying group-mean centering for listenBPM and responseBPM (by subject only)...")
df_intero <- df_intero %>%
  group_by(subject) %>%
  mutate(
    listenBPM_centered = listenBPM - mean(listenBPM, na.rm = TRUE),
    responseBPM_centered = responseBPM - mean(responseBPM, na.rm = TRUE)
  ) %>%
  ungroup()

# --- Model Fitting ---
# We are using RAW (but centered) variables.

message("Fitting model: Response ~ drugs * listenBPM_centered + responseBPM_centered + nTrials + visit + (1 + listenBPM_centered + responseBPM_centered | subject)")

# Helper for convergence check
is_converged <- function(m) {
  if (inherits(m, "try-error")) return(FALSE)
  if (m$fit$convergence != 0) return(FALSE)
  if (!m$sdr$pdHess) return(FALSE)
  return(TRUE)
}

mod_formula <- Response_binary ~ drugs * listenBPM_centered + responseBPM_centered + nTrials + visit

# 1. Maximal Random Effects
mod <- try(glmmTMB(
  formula = update(mod_formula, . ~ . + (1 + listenBPM_centered + responseBPM_centered | subject)),
  data = df_intero,
  family = binomial(link = "logit"),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
), silent = TRUE)

# 2. Uncorrelated Random Effects
if (!is_converged(mod)) {
  message("Maximal model failed. Trying uncorrelated random effects (||)...")
  mod <- try(glmmTMB(
    formula = update(mod_formula, . ~ . + (1 + listenBPM_centered + responseBPM_centered || subject)),
    data = df_intero,
    family = binomial(link = "logit")
  ), silent = TRUE)
}

# 3. Random Intercepts Only
if (!is_converged(mod)) {
  message("Random slopes failed. Dropping random slopes...")
  mod <- try(glmmTMB(
    formula = update(mod_formula, . ~ . + (1 | subject)),
    data = df_intero,
    family = binomial(link = "logit")
  ), silent = TRUE)
}

if (!is_converged(mod)) {
  message("WARNING: Model failed to converge.")
} else {
  message("Model converged successfully.")
}

# --- Summary & Saving (concise) ---
mod_sum <- summary(mod)

# Save model (RDS)
model_filename <- sprintf("model_raw_interaction_%s.rds", TIMESTAMP)
model_file <- file.path(DATA_DIR, model_filename)
saveRDS(mod, model_file)
message("Saved model to: ", model_file)

# Save fixed-effects coefficients as CSV (concise)
coef_df <- tryCatch({
  as.data.frame(mod_sum$coefficients$cond) %>% tibble::rownames_to_column(var = "term")
}, error = function(e) {
  message("Could not extract coefficients table: ", e$message)
  NULL
})
if (!is.null(coef_df)) {
  coef_csv <- file.path(DOCS_DIR, sprintf("model_coefficients_%s.csv", TIMESTAMP))
  readr::write_csv(coef_df, coef_csv)
  message("Saved coefficients to: ", coef_csv)
  message("Top coefficients:")
  print(utils::head(coef_df, 10))
}

# Save full summary to text using capture.output (avoids flooding console)
summary_filename <- sprintf("model_summary_raw_%s.txt", TIMESTAMP)
summary_file <- file.path(DOCS_DIR, summary_filename)
tryCatch({
  capture.output(print(mod_sum), file = summary_file)
  message("Saved model summary to: ", summary_file)
}, error = function(e) {
  message("Failed to save full summary: ", e$message)
})

# --- Post-hoc Tests (safe, saved to files) ---
message("Calculating simple slopes for listenBPM by Drug (emm_trends)...")
emm_trends <- tryCatch({
  emmeans::emtrends(mod, specs = "drugs", var = "listenBPM_centered")
}, error = function(e) {
  message("Emtrends failed: ", e$message)
  NULL
})
if (!is.null(emm_trends)) {
  emm_trends_df <- as.data.frame(summary(emm_trends))
  emm_filename <- file.path(DATA_DIR, sprintf("emmeans_listenBPM_slopes_raw_%s.csv", TIMESTAMP))
  readr::write_csv(emm_trends_df, emm_filename)
  message("Saved emmeans trends to: ", emm_filename)
  message("Emmeans trends (top):")
  print(utils::head(emm_trends_df, 10))

  # Pairwise comparisons (safe)
  emm_pairs <- tryCatch({
    pairs(emm_trends)
  }, error = function(e) {
    message("Pairwise comparisons failed: ", e$message)
    NULL
  })
  if (!is.null(emm_pairs)) {
    pairs_df <- as.data.frame(summary(emm_pairs))
    pairs_filename <- file.path(DATA_DIR, sprintf("emmeans_pairs_listenBPM_%s.csv", TIMESTAMP))
    readr::write_csv(pairs_df, pairs_filename)
    message("Saved emmeans pairwise comparisons to: ", pairs_filename)
    print(utils::head(pairs_df, 10))
  }
}

message("Analysis complete!")
