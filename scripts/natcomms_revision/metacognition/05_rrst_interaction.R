#!/usr/bin/env Rscript
# Fit RRST confidence model with Stimulus * Drug * Accuracy three-way interaction
# Examines the effect of allowing stimulus intensity to covary with response accuracy
# Stimulus is subject-mean-centered before entering the model
# Model: Conf ~ Stimulus_c * Drug * Accuracy + (1 + Stimulus_c + Accuracy | subject)

if (!require("pacman")) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(dplyr, readr, tidyr, glmmTMB, emmeans, ggplot2, ggeffects, sjPlot, here)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision", "metacognition")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
DOCS_DIR <- file.path(OUTPUT_DIR, "docs")

FIGS_DIR <- file.path(OUTPUT_DIR, "figs")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Run Timestamp: ", TIMESTAMP)

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

# --- Load and preprocess RRST data (identical to 01_model.R) ---
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

RRST_trial_data <- clean_and_report(RRST_trial_data, "RRST_trial_data")

# Subject-level mean-centering of Stimulus
RRST_trial_data <- RRST_trial_data %>%
  group_by(subject) %>%
  mutate(Stimulus_c = Stimulus - mean(Stimulus, na.rm = TRUE)) %>%
  ungroup()

message("RRST subjects (before exclusion): ", length(unique(RRST_trial_data$subject)))
message("RRST trials  (before exclusion): ", nrow(RRST_trial_data))

# --- Exclude random-effect outlier subjects ---
# Identified via 06_rrst_interaction_diagnostics.R: subjects with random effects
# > 2 SD from the group mean on one or more terms.
#   sub_3004  - outlier on Stimulus_c slope (steep positive)
#   sub_3005  - outlier on Accuracy intercept (large negative)
#   sub_4032  - outlier on Accuracy intercept (large positive)
#   sub_4054  - outlier on Intercept (very high baseline confidence)
#   sub_4068  - outlier on 3 terms: Stimulus_c, Accuracy, Stimulus_c:Accuracy
#   sub_4076  - outlier on Accuracy intercept (large negative)
#   sub_4079  - outlier on Intercept and Stimulus_c slope
RE_OUTLIER_SUBJECTS <- c("sub_3004", "sub_3005", "sub_4032", "sub_4054",
                         "sub_4068", "sub_4076", "sub_4079")

RRST_trial_data <- RRST_trial_data %>%
  filter(!subject %in% RE_OUTLIER_SUBJECTS)

message("Excluded ", length(RE_OUTLIER_SUBJECTS), " RE-outlier subjects: ",
        paste(RE_OUTLIER_SUBJECTS, collapse = ", "))

# --- Filter to Stimulus >= 8 (drop sparse low levels 5-7) ---
RRST_trial_data <- RRST_trial_data %>%
  filter(Stimulus >= 8)

message("Filtered to Stimulus >= 8")
message("RRST subjects (after exclusions): ", length(unique(RRST_trial_data$subject)))
message("RRST trials  (after exclusions): ", nrow(RRST_trial_data))

# --- Fit three-way interaction model ---
message("Fitting RRST Stimulus_c * Drug * Accuracy model...")

rrst_conf_model_inter <- glmmTMB(
  Conf ~ Stimulus_c * Drug * Accuracy + (1 + Stimulus_c * Accuracy || subject),
  data = RRST_trial_data,
  family = ordbeta(),
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)

mod_sum <- summary(rrst_conf_model_inter)

message("\nFixed effects:")
print(mod_sum$coefficients$cond)
message("\nAIC: ", AIC(rrst_conf_model_inter))

# Save model
model_filename <- sprintf("model_rrst_confidence_interaction_%s.rds", TIMESTAMP)
model_file <- file.path(DATA_DIR, model_filename)
saveRDS(rrst_conf_model_inter, model_file)
message("Saved model to: ", model_file)

# Save summary
summary_filename <- sprintf("model_summary_rrst_confidence_interaction_%s.txt", TIMESTAMP)
summary_file <- file.path(DOCS_DIR, summary_filename)
sink(summary_file)
cat("========================================\n")
cat("Model: RRST Confidence - Stimulus * Drug * Accuracy Interaction\n")
cat("Run Timestamp: ", TIMESTAMP, "\n")
cat("Formula: Conf ~ Stimulus_c * Drug * Accuracy + (1 + Stimulus_c * Accuracy || subject)\n")
cat("Note: Stimulus_c is subject-mean-centered Stimulus intensity\n")
cat("Family: ordbeta()\n")
cat("Subjects: ", length(unique(RRST_trial_data$subject)), "\n")
cat("Trials: ", nrow(RRST_trial_data), "\n")
cat("RE-outlier subjects excluded: ", paste(RE_OUTLIER_SUBJECTS, collapse = ", "), "\n")
cat("========================================\n\n")
print(mod_sum)
cat("\n\nAIC: ", AIC(rrst_conf_model_inter), "\n")
sink()
message("Saved summary to: ", summary_file)

# ============================================================
# Plots
# ============================================================
message("\n--- Generating plots ---")

# --- 1. Stimulus slopes by Drug x Accuracy (emtrends) ---
# Shows how the Stimulus-confidence slope differs across Drug and Accuracy
message("Computing Stimulus slopes by Drug x Accuracy...")
stim_trends <- emtrends(rrst_conf_model_inter, ~ Drug | Accuracy,
                        var = "Stimulus_c")
message("Stimulus slopes (emtrends):")
print(summary(stim_trends))

# Pairwise contrasts of slopes within each Accuracy level
stim_trend_pairs <- pairs(stim_trends)
message("Pairwise slope contrasts:")
print(summary(stim_trend_pairs))

# Plot the slopes
trend_df <- as.data.frame(summary(stim_trends))
z_val <- qnorm(0.975)
trend_df$lower <- trend_df$Stimulus_c.trend - z_val * trend_df$SE
trend_df$upper <- trend_df$Stimulus_c.trend + z_val * trend_df$SE

p_slopes <- ggplot(trend_df,
                   aes(x = Drug, y = Stimulus_c.trend,
                       color = Drug)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8,
                  linewidth = 1.0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ Accuracy) +
  labs(y = "Stimulus Slope (emtrend)",
       x = "",
       title = "Stimulus-Confidence Slopes by Drug and Accuracy") +
  scale_color_manual(values = c("placebo" = "#2E3FAE",
                                "propranolol" = "#C13A8C",
                                "bisoprolol" = "#1FA187")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 13))

ggsave(file.path(FIGS_DIR,
                 sprintf("rrst_interaction_emtrends_%s.png", TIMESTAMP)),
       p_slopes, width = 8, height = 5, dpi = 300, bg = "white")
message("Saved emtrends plot.")

# --- 2. Predicted confidence across Stimulus by Drug, faceted by Accuracy ---
# Uses sjPlot::plot_model for the continuous x categorical interaction
p_pred <- plot_model(rrst_conf_model_inter, type = "pred",
                     terms = c("Stimulus_c [all]", "Drug", "Accuracy"),
                     colors = c("#2E3FAE", "#C13A8C", "#1FA187")) +
  labs(x = "Stimulus Intensity (subject-centered)",
       y = "Predicted Confidence",
       title = "Predicted Confidence: Stimulus x Drug x Accuracy",
       color = "Drug") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 13))

ggsave(file.path(FIGS_DIR,
                 sprintf("rrst_interaction_predicted_%s.png", TIMESTAMP)),
       p_pred, width = 10, height = 5, dpi = 300, bg = "white")
message("Saved predicted confidence plot.")

# --- 3. Emmeans at low/medium/high Stimulus by Drug x Accuracy ---
# Evaluates predicted confidence at -1 SD, 0, +1 SD of Stimulus_c
stim_sd <- sd(RRST_trial_data$Stimulus_c)
stim_levels <- c(-stim_sd, 0, stim_sd)
stim_labels <- c("Low (-1 SD)", "Mean", "High (+1 SD)")

emm_stim <- emmeans(rrst_conf_model_inter,
                    ~ Drug * Accuracy | Stimulus_c,
                    at = list(Stimulus_c = stim_levels))
emm_df <- as.data.frame(summary(emm_stim))
emm_df$Stimulus_level <- factor(
  rep(stim_labels, each = nrow(emm_df) / 3),
  levels = stim_labels
)

p_emm <- ggplot(emm_df,
                aes(x = Drug, y = emmean, color = Accuracy,
                    group = Accuracy)) +
  geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                  position = position_dodge(width = 0.4),
                  size = 0.7, linewidth = 0.9) +
  facet_wrap(~ Stimulus_level) +
  scale_color_manual(values = c("Correct" = "#18699D",
                                "Incorrect" = "#E8BB1B")) +
  labs(y = "Estimated Confidence",
       x = "",
       title = "Drug x Accuracy at Low / Mean / High Stimulus") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 13))

ggsave(file.path(FIGS_DIR,
                 sprintf("rrst_interaction_emmeans_stimulus_%s.png", TIMESTAMP)),
       p_emm, width = 10, height = 5, dpi = 300, bg = "white")
message("Saved emmeans at stimulus levels plot.")

# --- 4. Model vs raw data: Accuracy ribbons faceted by Drug ---
# Adapted from original Analysis.Rmd plot_3 style
message("Generating model vs data plot (Accuracy x Drug)...")

P1 <- ggeffects::ggemmeans(rrst_conf_model_inter,
                           terms = c("Stimulus_c", "Accuracy", "Drug"))

# Group by original discrete Stimulus levels (not continuous Stimulus_c)
# then compute mean Stimulus_c per group for x-axis alignment with model
subj_data <- RRST_trial_data %>%
  group_by(Accuracy, Stimulus, Drug) %>%
  summarise(
    mean = mean(Conf),
    q5  = mean(Conf) - 2 * sd(Conf) / sqrt(n()),
    q95 = mean(Conf) + 2 * sd(Conf) / sqrt(n()),
    Stimulus_c = mean(Stimulus_c),
    .groups = "drop"
  ) %>%
  mutate(facet = Drug)

p_model_data <- data.frame(P1) %>%
  mutate(Accuracy = as.factor(group)) %>%
  ggplot() +
  geom_line(aes(x = x, y = predicted, color = Accuracy)) +
  geom_ribbon(aes(x = x, y = predicted,
                  ymin = conf.low, ymax = conf.high,
                  fill = Accuracy),
              alpha = 0.25) +
  geom_pointrange(data = subj_data,
                  aes(x = Stimulus_c, y = mean,
                      ymin = q5, ymax = q95,
                      fill = as.factor(Accuracy)),
                  shape = 21, col = "black") +
  scale_color_manual(values = c("Correct" = "green", "Incorrect" = "red")) +
  scale_fill_manual(values = c("Correct" = "green", "Incorrect" = "red")) +
  facet_wrap(~ facet) +
  labs(x = "Stimulus Intensity (subject-centered)",
       y = "Confidence",
       title = "Model Predictions vs Raw Data by Drug") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 13))

ggsave(file.path(FIGS_DIR,
                 sprintf("rrst_interaction_model_vs_data_%s.png", TIMESTAMP)),
       p_model_data, width = 12, height = 5, dpi = 300, bg = "white")
message("Saved model vs data plot.")

message("\nDone!")
