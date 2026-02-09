#!/usr/bin/env Rscript
# Plot model predictions vs raw data for the RRST interaction model

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, tidyr, glmmTMB, ggplot2, ggeffects, here)

BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision", "metacognition")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

# --- Load model ---
model_files <- list.files(DATA_DIR,
                          pattern = "model_rrst_confidence_interaction_.*\\.rds",
                          full.names = TRUE)
details <- file.info(model_files)
model_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
message("Loading model: ", basename(model_file))
mod <- readRDS(model_file)

# --- Recreate data (same preprocessing + exclusions as 05_rrst_interaction.R) ---
RE_OUTLIER_SUBJECTS <- c("sub_3004", "sub_3005", "sub_4032", "sub_4054",
                         "sub_4068", "sub_4076", "sub_4079")

RRST_trial_data <- readr::read_csv(
  file.path(BASE_DIR, "data", "cleaned", "RRST.csv"),
  show_col_types = FALSE
) %>%
  filter(ConfResp >= 0 & ConfResp <= 100) %>%
  mutate(
    Conf = ConfResp / 100,
    Stimulus = Stim,
    Drug = factor(drugs,
                  levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    Accuracy = factor(Resp, levels = c(1, 0), labels = c("Correct", "Incorrect"))
  ) %>%
  drop_na() %>%
  group_by(subject) %>%
  mutate(Stimulus_c = Stimulus - mean(Stimulus, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!subject %in% RE_OUTLIER_SUBJECTS)

message("Data: ", nrow(RRST_trial_data), " trials, ",
        length(unique(RRST_trial_data$subject)), " subjects")

# --- Summarise raw data by discrete Stimulus level ---
# Stimulus has 13 natural levels (5-17), so group by those
# then use mean Stimulus_c per group for x-axis position
raw_summary <- RRST_trial_data %>%
  group_by(Drug, Accuracy, Stimulus) %>%
  summarise(
    Conf_mean = mean(Conf, na.rm = TRUE),
    Conf_se = sd(Conf, na.rm = TRUE) / sqrt(n()),
    Stimulus_c = mean(Stimulus_c, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

message("Raw data summary: ", nrow(raw_summary), " points (",
        length(unique(raw_summary$Stimulus)), " stimulus levels x ",
        "3 drugs x 2 accuracy)")

# --- Get model predictions ---
preds <- ggpredict(mod, terms = c("Stimulus_c [n=100]", "Drug", "Accuracy"))

preds <- as.data.frame(preds)
names(preds)[names(preds) == "group"] <- "Drug"
names(preds)[names(preds) == "facet"] <- "Accuracy"

# --- Drug colours ---
drug_cols <- c("placebo" = "#2E3FAE",
               "propranolol" = "#C13A8C",
               "bisoprolol" = "#1FA187")

# --- Plot: faceted by Accuracy, colored by Drug ---
p <- ggplot() +
  geom_ribbon(
    data = preds,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = Drug),
    alpha = 0.15
  ) +
  geom_line(
    data = preds,
    aes(x = x, y = predicted, color = Drug),
    linewidth = 1.0
  ) +
  geom_pointrange(
    data = raw_summary,
    aes(x = Stimulus_c, y = Conf_mean,
        ymin = Conf_mean - Conf_se, ymax = Conf_mean + Conf_se,
        color = Drug),
    size = 0.3, linewidth = 0.4, alpha = 0.7
  ) +
  facet_wrap(~ Accuracy, scales = "fixed") +
  scale_color_manual(values = drug_cols) +
  scale_fill_manual(values = drug_cols) +
  labs(
    x = "Stimulus Intensity (subject-centered)",
    y = "Confidence",
    color = "Drug", fill = "Drug",
    title = "Model Predictions vs Raw Data (Binned by Stimulus Level)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 14)
  )

outfile <- file.path(FIGS_DIR,
                     sprintf("rrst_interaction_model_vs_data_%s.png", TIMESTAMP))
ggsave(outfile, p, width = 10, height = 5, dpi = 300, bg = "white")
message("Saved to: ", outfile)
