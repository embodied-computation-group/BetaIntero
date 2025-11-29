#!/usr/bin/env Rscript
# Fit hierarchical model predicting choice from deltaBPM (listenBPM - responseBPM)
# Model: Decision ~ drugs * deltaBPM
# Random Effects: Maximal random slopes for deltaBPM
# Excludes outlier sub_4071

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
  mutate(drugs = factor(drugs, levels = c("PLACEBO", "BISO", "PROP")))

message("Interoceptive trials with complete data: ", nrow(df_intero))

# --- Outlier Removal ---
outliers <- c("sub_4071")
message("Removing known outlier subjects: ", paste(outliers, collapse = ", "))
df_intero <- df_intero %>%
  filter(!subject %in% outliers)

# --- Feature Engineering: Diff BPM ---
# Task: "Are tones (response) faster than heart (listen)?"
# Predictor: diffBPM = responseBPM - listenBPM
# If diffBPM > 0 (Response > Listen), Tones are Faster -> Expect "More" (1)
# We expect a POSITIVE slope.

df_intero <- df_intero %>%
  mutate(diffBPM = responseBPM - listenBPM)

# Ensure visit is a factor
df_intero <- df_intero %>%
  mutate(visit = as.factor(visit))

# --- Scaling ---
# Scale nTrials within subject-drug session
safe_scale <- function(x) {
  if(sd(x, na.rm=TRUE) == 0) return(rep(0, length(x)))
  return(as.vector(scale(x)))
}

df_intero <- df_intero %>%
  group_by(subject, drugs) %>%
  mutate(
    nTrials_scaled = safe_scale(nTrials)
  ) %>%
  ungroup()

# --- Model Fitting ---
message("Fitting model: Response ~ drugs * diffBPM + nTrials_scaled + visit + (1 + diffBPM | subject)")

# Helper for convergence check
is_converged <- function(m) {
  if (inherits(m, "try-error")) return(FALSE)
  if (m$fit$convergence != 0) return(FALSE)
  if (!m$sdr$pdHess) return(FALSE)
  return(TRUE)
}

# 1. Maximal Random Effects (Unscaled diffBPM)
mod <- try(glmmTMB(
  Response_binary ~ drugs * diffBPM + nTrials_scaled + visit +
    (1 + diffBPM | subject),
  data = df_intero,
  family = binomial(link = "logit")
), silent = TRUE)

# 2. Uncorrelated Random Effects
if (!is_converged(mod)) {
  message("Maximal model failed. Trying uncorrelated random effects (||)...")
  mod <- try(glmmTMB(
    Response_binary ~ drugs * diffBPM + nTrials_scaled + visit +
      (1 + diffBPM || subject),
    data = df_intero,
    family = binomial(link = "logit")
  ), silent = TRUE)
}

# 3. Random Intercepts Only
if (!is_converged(mod)) {
  message("Random slopes failed. Dropping random slopes...")
  mod <- try(glmmTMB(
    Response_binary ~ drugs * diffBPM + nTrials_scaled + visit +
      (1 | subject),
    data = df_intero,
    family = binomial(link = "logit")
  ), silent = TRUE)
}

if (!is_converged(mod)) {
  message("WARNING: Model failed to converge.")
} else {
  message("Model converged successfully.")
}

# --- Summary & Saving ---
mod_sum <- summary(mod)
print(mod_sum$coefficients$cond)

# Save model
model_filename <- sprintf("model_diffBPM_interaction_%s.rds", TIMESTAMP)
model_file <- file.path(DATA_DIR, model_filename)
saveRDS(mod, model_file)
message("Saved model to: ", model_file)

# Save summary text
summary_filename <- sprintf("model_summary_diffBPM_%s.txt", TIMESTAMP)
summary_file <- file.path(DOCS_DIR, summary_filename)
sink(summary_file)
print(mod_sum)
sink()

# --- Calculate Thresholds (PSE) ---
message("\nCalculating Thresholds (PSE) for each drug...")

# Extract fixed effects
fixefs <- fixef(mod)$cond
vcovs <- vcov(mod)$cond

# Helper to get PSE and CI via simulation
get_pse <- function(drug_name, n_sim = 10000) {
  # Simulate coefficients
  sim_coefs <- MASS::mvrnorm(n_sim, fixefs, vcovs)
  
  if (drug_name == "PLACEBO") {
    b0 <- sim_coefs[, "(Intercept)"]
    b1 <- sim_coefs[, "diffBPM"]
  } else {
    b0 <- sim_coefs[, "(Intercept)"] + sim_coefs[, paste0("drugs", drug_name)]
    b1 <- sim_coefs[, "diffBPM"] + sim_coefs[, paste0("drugs", drug_name, ":diffBPM")]
  }
  
  # PSE = -b0 / b1
  pses <- -b0 / b1
  
  return(c(
    mean = mean(pses),
    lower = quantile(pses, 0.025),
    upper = quantile(pses, 0.975)
  ))
}

pse_res <- data.frame(
  drugs = c("PLACEBO", "BISO", "PROP"),
  PSE = NA,
  CI_Lower = NA,
  CI_Upper = NA
)

pse_res[1, 2:4] <- get_pse("PLACEBO")
pse_res[2, 2:4] <- get_pse("BISO")
pse_res[3, 2:4] <- get_pse("PROP")

print(pse_res)

# Save PSE results
pse_filename <- sprintf("pse_thresholds_diffBPM_%s.csv", TIMESTAMP)
readr::write_csv(pse_res, file.path(DATA_DIR, pse_filename))
message("Saved PSE thresholds to: ", file.path(DATA_DIR, pse_filename))

# --- Post-hoc Tests ---
message("\nCalculating simple slopes for diffBPM by Drug...")
emm_trends <- emmeans::emtrends(mod, specs = "drugs", var = "diffBPM")
print(summary(emm_trends))

# Save Emmeans
emm_filename <- sprintf("emmeans_diffBPM_slopes_%s.csv", TIMESTAMP)
readr::write_csv(as.data.frame(summary(emm_trends)), file.path(DATA_DIR, emm_filename))

message("Analysis complete!")
