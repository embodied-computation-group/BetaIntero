#!/usr/bin/env Rscript
# Fit simple model: Response ~ listenBPM * visit + responseBPM

if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org")
pacman::p_load(dplyr, readr, glmmTMB, emmeans, here)

# Directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "revision1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
DOCS_DIR <- file.path(OUTPUT_DIR, "docs")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Run Timestamp: ", TIMESTAMP)

# Load data
INPUT_CSV <- file.path(BASE_DIR, "data", "HRD_trial_level_data.csv")
if (!file.exists(INPUT_CSV)) stop("Input CSV not found: ", INPUT_CSV)

message("Reading data from: ", INPUT_CSV)
df <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)

# Filter & clean
if ("Modality" %in% names(df)) {
  df_intero <- df %>% filter(Modality == "Intero")
} else {
  df_intero <- df
}

df_intero <- df_intero %>%
  filter(Decision %in% c("More", "Less")) %>%
  filter(!is.na(listenBPM), !is.na(responseBPM)) %>%
  mutate(Response_binary = ifelse(Decision == "More", 1, 0)) %>%
  filter(drugs %in% c("PLACEBO", "BISO", "PROP")) %>%
  mutate(drugs = factor(drugs, levels = c("PLACEBO", "BISO", "PROP"))) %>%
  mutate(visit = as.factor(visit))

# Remove known outliers
outliers <- c("sub_4071")
df_intero <- df_intero %>% filter(!subject %in% outliers)

message("Observations after filtering: ", nrow(df_intero))

# Center within-subject (subject-only centering)
df_intero <- df_intero %>%
  group_by(subject) %>%
  mutate(
    listenBPM_centered = listenBPM - mean(listenBPM, na.rm = TRUE),
    responseBPM_centered = responseBPM - mean(responseBPM, na.rm = TRUE)
  ) %>%
  ungroup()

# Model formula
mod_formula <- Response_binary ~ listenBPM_centered * visit + responseBPM_centered
message("Fitting model: Response ~ listenBPM_centered * visit + responseBPM_centered + (random effects)")

# Convergence helper
is_converged <- function(m) {
  if (inherits(m, "try-error")) return(FALSE)
  if (!is.null(m$fit) && !is.null(m$fit$convergence) && m$fit$convergence != 0) return(FALSE)
  if (!is.null(m$sdr) && !is.null(m$sdr$pdHess) && !m$sdr$pdHess) return(FALSE)
  return(TRUE)
}

# 1. Maximal random effects
mod <- try(glmmTMB(
  formula = update(mod_formula, . ~ . + (1 + listenBPM_centered + responseBPM_centered | subject)),
  data = df_intero,
  family = binomial(link = "logit"),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
), silent = TRUE)

# 2. Uncorrelated random slopes
if (!is_converged(mod)) {
  message("Maximal model failed; trying uncorrelated random slopes...")
  mod <- try(glmmTMB(
    formula = update(mod_formula, . ~ . + (1 + listenBPM_centered + responseBPM_centered || subject)),
    data = df_intero,
    family = binomial(link = "logit")
  ), silent = TRUE)
}

# 3. Random intercepts only
if (!is_converged(mod)) {
  message("Uncorrelated slopes failed; trying random intercepts only...")
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

# Save model
model_filename <- sprintf("model_listen_visit_%s.rds", TIMESTAMP)
model_file <- file.path(DATA_DIR, model_filename)
saveRDS(mod, model_file)
message("Saved model to: ", model_file)

# Save concise coefficients
mod_sum <- summary(mod)
coef_df <- tryCatch({
  as.data.frame(mod_sum$coefficients$cond) %>% tibble::rownames_to_column(var = "term")
}, error = function(e) {
  message("Failed to extract coefficients: ", e$message)
  NULL
})
if (!is.null(coef_df)) {
  coef_csv <- file.path(DOCS_DIR, sprintf("model_listen_visit_coefficients_%s.csv", TIMESTAMP))
  readr::write_csv(coef_df, coef_csv)
  message("Saved coefficients to: ", coef_csv)
  print(utils::head(coef_df, 10))
}

# Save summary
summary_file <- file.path(DOCS_DIR, sprintf("model_listen_visit_summary_%s.txt", TIMESTAMP))
tryCatch({
  capture.output(print(mod_sum), file = summary_file)
  message("Saved model summary to: ", summary_file)
}, error = function(e) {
  message("Failed to save summary: ", e$message)
})

# Post-hoc: trends of listenBPM by visit
emm_trends <- tryCatch({
  emmeans::emtrends(mod, specs = "visit", var = "listenBPM_centered")
}, error = function(e) {
  message("Emtrends failed: ", e$message)
  NULL
})
if (!is.null(emm_trends)) {
  emm_df <- as.data.frame(summary(emm_trends))
  emm_file <- file.path(DATA_DIR, sprintf("emmeans_listenBPM_by_visit_%s.csv", TIMESTAMP))
  readr::write_csv(emm_df, emm_file)
  message("Saved emmeans trends to: ", emm_file)
}

message("Done.")
