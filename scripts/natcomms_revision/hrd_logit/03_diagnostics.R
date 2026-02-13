#!/usr/bin/env Rscript
# Diagnostics for the Interoception Interaction Model
# 1. Check raw data patterns
# 2. Check collinearity
# 3. Check for influential subjects

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, ggplot2, glmmTMB, here, tidyr, patchwork, performance)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision", "hrd_logit")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Diagnostics Run Timestamp: ", TIMESTAMP)

# Load data
INPUT_CSV <- file.path(BASE_DIR, "data", "HRD_trial_level_data.csv")
df <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)

# Preprocessing (same as model script)
cleaned_subjects <- unique(read.csv(file.path(BASE_DIR, "data", "cleaned", "HRD.csv"))$subject)

df_intero <- df %>%
  filter(Modality == "Intero") %>%
  filter(Decision %in% c("More", "Less")) %>%
  filter(!is.na(listenBPM), !is.na(responseBPM)) %>%
  mutate(Response_binary = ifelse(Decision == "More", 1, 0)) %>%
  filter(drugs %in% c("PLACEBO", "BISO", "PROP")) %>%
  mutate(drugs = factor(drugs, levels = c("PLACEBO", "BISO", "PROP"))) %>%
  filter(subject %in% cleaned_subjects)

message("Data loaded. N = ", nrow(df_intero), " (", length(unique(df_intero$subject)), " subjects)")

# --- 1. Check Correlations (Collinearity) ---
message("\n--- Checking Correlations ---")
cor_df <- df_intero %>%
  group_by(drugs) %>%
  summarise(
    cor_listen_response = cor(listenBPM, responseBPM),
    cor_listen_choice = cor(listenBPM, Response_binary),
    cor_response_choice = cor(responseBPM, Response_binary)
  )
print(cor_df)

# Plot correlation between listenBPM and responseBPM by drug
p_cor <- ggplot(df_intero, aes(x = listenBPM, y = responseBPM, color = drugs)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm") +
  facet_wrap(~drugs) +
  theme_minimal() +
  labs(title = "Correlation: Listen BPM vs Response BPM")

ggsave(file.path(FIGS_DIR, sprintf("diag_correlation_listen_response_%s.png", TIMESTAMP)), p_cor, width = 10, height = 4, bg = "white")


# --- 2. Visualizing the Raw Effect ---
message("\n--- Visualizing Raw Effect ---")
# Bin listenBPM and calculate probability of "More"
df_binned <- df_intero %>%
  mutate(listen_bin = cut(listenBPM, breaks = 10)) %>%
  group_by(drugs, listen_bin) %>%
  summarise(
    prob_more = mean(Response_binary),
    n = n(),
    listen_mean = mean(listenBPM),
    .groups = "drop"
  )

p_raw <- ggplot(df_binned, aes(x = listen_mean, y = prob_more, color = drugs, group = drugs)) +
  geom_point(aes(size = n)) +
  geom_line() +
  geom_smooth(data = df_intero, aes(x = listenBPM, y = Response_binary), 
              method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  theme_minimal() +
  labs(title = "Raw Data: Probability of 'More' vs Listen BPM",
       y = "Pr(Response = More)", x = "Listen BPM")

ggsave(file.path(FIGS_DIR, sprintf("diag_raw_effect_%s.png", TIMESTAMP)), p_raw, width = 8, height = 6, bg = "white")

# --- 2b. Raw Effect Controlling for Response BPM ---
message("\n--- Visualizing Raw Effect (controlling responseBPM) ---")
# Residualize listenBPM: regress out responseBPM within each drug condition
# so the x-axis reflects listenBPM variance independent of responseBPM
df_intero <- df_intero %>%
  group_by(drugs) %>%
  mutate(listenBPM_resid = residuals(lm(listenBPM ~ responseBPM))) %>%
  ungroup()

df_binned_adj <- df_intero %>%
  mutate(listen_bin = cut(listenBPM_resid, breaks = 10)) %>%
  group_by(drugs, listen_bin) %>%
  summarise(
    prob_more = mean(Response_binary),
    n = n(),
    listen_mean = mean(listenBPM_resid),
    .groups = "drop"
  )

p_raw_adj <- ggplot(df_binned_adj, aes(x = listen_mean, y = prob_more, color = drugs, group = drugs)) +
  geom_point(aes(size = n)) +
  geom_line() +
  geom_smooth(data = df_intero, aes(x = listenBPM_resid, y = Response_binary),
              method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  theme_minimal() +
  labs(title = "Raw Data: Pr('More') vs Listen BPM (residualized on Response BPM)",
       subtitle = "listenBPM after removing variance shared with responseBPM",
       y = "Pr(Response = More)", x = "Listen BPM (residualized)")

ggsave(file.path(FIGS_DIR, sprintf("diag_raw_effect_adj_%s.png", TIMESTAMP)), p_raw_adj, width = 8, height = 6, bg = "white")
message("Saved adjusted raw effect plot.")

# --- 3. Influential Subjects Check ---
message("\n--- Checking Influential Subjects ---")
# We will calculate the slope of listenBPM -> Response for each subject individually using simple glm
# This helps identify if a few subjects are driving the drug difference

subj_slopes <- df_intero %>%
  group_by(subject, drugs) %>%
  do({
    mod <- glm(Response_binary ~ scale(listenBPM) + scale(responseBPM), data = ., family = binomial)
    data.frame(
      intercept = coef(mod)[1],
      slope_listen = coef(mod)[2],
      slope_response = coef(mod)[3],
      n_trials = nrow(.)
    )
  }) %>%
  ungroup()

# Plot distribution of slopes by drug
p_slopes <- ggplot(subj_slopes, aes(x = drugs, y = slope_listen, fill = drugs)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Distribution of Individual Subject Slopes (Listen BPM)",
       y = "Slope (Log Odds)")

ggsave(file.path(FIGS_DIR, sprintf("diag_subject_slopes_%s.png", TIMESTAMP)), p_slopes, width = 6, height = 6, bg = "white")

# Identify potential outliers (e.g., slope > 2 SD from mean)
slope_stats <- subj_slopes %>%
  group_by(drugs) %>%
  mutate(
    mean_slope = mean(slope_listen, na.rm = TRUE),
    sd_slope = sd(slope_listen, na.rm = TRUE),
    is_outlier = abs(slope_listen - mean_slope) > 2 * sd_slope
  )

outliers <- slope_stats %>% filter(is_outlier)
print("Potential outlier subjects based on individual slopes:")
print(outliers %>% select(subject, drugs, slope_listen))

# --- 4. Check VIF in the main model ---
message("\n--- Checking VIF ---")
# Find most recent model file matching the pattern
model_pattern <- "model_combined_drug_interaction_with_trialN_.*\\.rds"
model_files <- list.files(DATA_DIR, pattern = model_pattern, full.names = TRUE)

if (length(model_files) == 0) {
  # Fallback
  model_file <- file.path(DATA_DIR, "model_combined_drug_interaction_with_trialN.rds")
} else {
  details <- file.info(model_files)
  model_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
}

if (file.exists(model_file)) {
  message("Loading model for VIF check: ", basename(model_file))
  mod <- readRDS(model_file)
  vif_res <- performance::check_collinearity(mod)
  print(vif_res)
} else {
  message("Model file not found, skipping VIF check.")
}

message("Diagnostics complete.")
