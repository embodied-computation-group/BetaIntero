#!/usr/bin/env Rscript
# Plot metacognitive confidence model results (RRST + HRD)
# Reproduces plot_3() from scripts/plots.R using patchwork layout

# Ensure required packages are installed
if (!require("pacman")) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(dplyr, readr, glmmTMB, ggplot2, ggeffects, emmeans, here, patchwork)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision", "metacognition")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")

# Generate Timestamp
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Plot Run Timestamp: ", TIMESTAMP)

# --- Load models ---
load_latest_model <- function(pattern, data_dir) {
  model_files <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  if (length(model_files) == 0) {
    stop("No model files found matching pattern: ", pattern, " in ", data_dir)
  }
  details <- file.info(model_files)
  model_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
  message("Loading model from: ", basename(model_file))
  readRDS(model_file)
}

rrst_model <- load_latest_model("model_rrst_confidence_.*\\.rds", DATA_DIR)
hrd_model  <- load_latest_model("model_hrd_confidence_.*\\.rds", DATA_DIR)

# --- Reconstruct data (needed for emmeans on RDS-loaded models) ---
RRST_trial_data <- readr::read_csv(
  file.path(BASE_DIR, "data", "cleaned", "RRST.csv"),
  show_col_types = FALSE
) %>%
  dplyr::filter(ConfResp >= 0 & ConfResp <= 100) %>%
  dplyr::mutate(
    Conf = ConfResp / 100,
    Stimulus = Stim,
    Drug = factor(drugs,
                  levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    Accuracy = factor(Resp, levels = c(1, 0),
                      labels = c("Correct", "Incorrect"))
  ) %>%
  tidyr::drop_na() %>%
  dplyr::group_by(subject) %>%
  dplyr::mutate(Stimulus_c = Stimulus - mean(Stimulus, na.rm = TRUE)) %>%
  dplyr::ungroup()

HRD_trl_data <- read.csv(file.path(BASE_DIR, "data", "cleaned", "HRD.csv")) %>%
  dplyr::mutate(
    Conf = Confidence / 100,
    drugs = as.factor(drugs),
    Drug = factor(drugs,
                  levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    ResponseCorrect = factor(ResponseCorrect,
                             levels = c("True", "False"),
                             labels = c("Correct", "Incorrect"))
  ) %>%
  dplyr::filter(subject != "sub_4049") %>%
  tidyr::drop_na() %>%
  dplyr::group_by(subject) %>%
  dplyr::mutate(BPM_scaled = as.numeric(scale(listenBPM))) %>%
  dplyr::ungroup()

# --- Colour palette and constants ---
ribbon_col <- "#B0B0B0"
binary_colourmap <- c("#18699D", "#E8BB1B")
contrast_colors <- c("#404040", "#009F73", "#B00089")
box_color <- "#F2F2F2"
line_width <- 1.3
fontsize <- 18

message("Generating plots...")

# ============================================================
# HRD: Predicted Confidence Plot (p1)
# ============================================================
preds_hrd <- ggpredict(hrd_model, terms = c("Drug", "ResponseCorrect"))

p1 <- ggplot(preds_hrd, aes(x = x, y = predicted, group = group, color = group)) +
  geom_line(size = 1.3, show.legend = FALSE) +
  geom_point(size = 4, show.legend = FALSE) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = ribbon_col),
              alpha = 0.15, outline.type = "both", show.legend = FALSE) +
  scale_color_manual(values = binary_colourmap) +
  scale_fill_manual(values = ribbon_col) +
  labs(
    x = "",
    y = "Predicted Confidence",
    color = "Response Correct",
    fill = ""
  ) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(0.43, 0.71)) +
  scale_x_discrete(expand = c(0.1, 0)) +
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7), labels = c("0.5", "0.6", "0.7")) +
  theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6), "pt") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = fontsize),
    legend.text = element_text(size = fontsize),
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = fontsize),
    axis.text.x = element_text(size = fontsize),
    axis.text.y = element_text(size = fontsize),
    plot.title = element_text(size = fontsize),
    axis.line.y = element_line(linewidth = 1),
    axis.line.x = element_line(linewidth = 1)
  )

# ============================================================
# HRD: Pairwise Contrast Forest Plot - Correct Trials (p2)
# ============================================================
emm_hrd <- emmeans(hrd_model, ~ Drug | ResponseCorrect, data = HRD_trl_data)
contrast_hrd <- contrast(emm_hrd, method = "pairwise")
contrast_hrd_df <- as.data.frame(contrast_hrd)

# Filter to Correct trials
contrast_hrd_correct <- subset(contrast_hrd_df, ResponseCorrect == "Correct")

# Calculate 95% CIs using z-distribution (df = Inf)
z_value <- qnorm(0.975)
contrast_hrd_correct$lower.CL <- contrast_hrd_correct$estimate - z_value * contrast_hrd_correct$SE
contrast_hrd_correct$upper.CL <- contrast_hrd_correct$estimate + z_value * contrast_hrd_correct$SE

p2 <- contrast_hrd_correct %>%
  mutate(contrast = factor(contrast, levels = c("propranolol - bisoprolol",
                                                 "placebo - bisoprolol",
                                                 "placebo - propranolol"))) %>%
  mutate(y_val = as.numeric(as.factor(contrast))) %>%
  ggplot(aes(col = contrast, show.legend = FALSE)) +
  geom_point(aes(y = y_val, x = estimate), size = 5, show.legend = FALSE) +
  geom_pointrange(aes(y = y_val, x = estimate, xmin = lower.CL, xmax = upper.CL),
                  linewidth = 1.5, show.legend = FALSE) +
  labs(
    title = "",
    y = "",
    x = "Difference in Predicted Confidence"
  ) +
  coord_cartesian(xlim = c(-0.325, 0), ylim = c(0.5, 3.4)) +
  geom_segment(aes(x = 0, xend = 0, y = -2.3, yend = 8.3),
               linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE) +
  theme_classic(base_size = 12) +
  scale_y_continuous(breaks = c(1, 2, 3), labels = c("prop - biso", "plac - biso", "plac - prop")) +
  theme(strip.background = element_rect(fill = box_color, color = "white"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_color_manual(values = contrast_colors) +
  theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6), "pt") +
  theme(
    plot.title = element_text(face = "bold", size = fontsize),
    axis.text.y = element_text(angle = 25, hjust = 0.5, size = fontsize),
    axis.text.x = element_text(size = fontsize),
    axis.title.x = element_text(size = fontsize),
    axis.line.x = element_line(linewidth = 1),
    legend.position = "none"
  )

# ============================================================
# RRST: Predicted Confidence Plot (p3)
# ============================================================
preds_rrst <- ggpredict(rrst_model, terms = c("Drug", "Accuracy"))
preds_rrst$group <- factor(preds_rrst$group, levels = c("Correct", "Incorrect"))

p3 <- ggplot(preds_rrst, aes(x = x, y = predicted, group = group, color = group)) +
  geom_line(size = 1.3) +
  geom_point(size = 4) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = ribbon_col),
              alpha = 0.15, outline.type = "both", show.legend = FALSE) +
  scale_color_manual(values = binary_colourmap) +
  scale_fill_manual(values = ribbon_col) +
  labs(
    x = "",
    y = "",
    color = "Accuracy",
    fill = ""
  ) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(0.43, 0.71)) +
  scale_x_discrete(expand = c(0.1, 0)) +
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7), labels = c("0.5", "0.6", "0.7")) +
  theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6), "pt") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = fontsize),
    legend.text = element_text(size = fontsize),
    axis.text = element_text(size = fontsize),
    axis.text.x = element_text(size = fontsize),
    axis.title = element_blank(),
    plot.title = element_text(size = 12, face = "bold"),
    axis.line.y = element_line(linewidth = 1),
    axis.line.x = element_line(linewidth = 1)
  )

# ============================================================
# RRST: Pairwise Contrast Forest Plot - Correct Trials (p4)
# ============================================================
emm_rrst <- emmeans(rrst_model, ~ Drug | Accuracy, data = RRST_trial_data)
contrast_rrst <- contrast(emm_rrst, method = "pairwise")
contrast_rrst_df <- as.data.frame(contrast_rrst)

# Filter to Correct trials
contrast_rrst_correct <- subset(contrast_rrst_df, Accuracy == "Correct")

# Calculate 95% CIs
contrast_rrst_correct$lower.CL <- contrast_rrst_correct$estimate - z_value * contrast_rrst_correct$SE
contrast_rrst_correct$upper.CL <- contrast_rrst_correct$estimate + z_value * contrast_rrst_correct$SE

p4 <- contrast_rrst_correct %>%
  mutate(contrast = factor(contrast, levels = c("propranolol - bisoprolol",
                                                 "placebo - bisoprolol",
                                                 "placebo - propranolol"))) %>%
  mutate(y_val = as.numeric(as.factor(contrast))) %>%
  ggplot(aes(col = contrast, show.legend = FALSE)) +
  geom_point(aes(y = y_val, x = estimate), size = 5, show.legend = FALSE) +
  geom_pointrange(aes(y = y_val, x = estimate, xmin = lower.CL, xmax = upper.CL),
                  linewidth = 1.5, show.legend = FALSE) +
  labs(
    title = "",
    y = "",
    x = "Difference in Predicted Confidence"
  ) +
  coord_cartesian(xlim = c(-0.1, 0.1), ylim = c(0.5, 3.4)) +
  geom_segment(aes(x = 0, xend = 0, y = -2.3, yend = 8.3),
               linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE) +
  theme_classic(base_size = 12) +
  scale_y_continuous(breaks = c(1, 2, 3), labels = c("", "", "")) +
  theme(strip.background = element_rect(fill = box_color, color = "white"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_color_manual(values = contrast_colors) +
  theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6), "pt") +
  theme(
    plot.title = element_text(face = "bold", size = fontsize),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = fontsize),
    axis.text.x = element_text(size = fontsize),
    axis.title.x = element_text(size = fontsize),
    axis.line.x = element_line(linewidth = 1),
    legend.position = "none"
  )

# ============================================================
# Combined Plot Assembly
# ============================================================
message("Assembling combined plot...")

# Top row: predicted confidence plots with shared legend
p_top <- (p1 + p3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Bottom row: pairwise contrast forest plots with title annotation
p_bottom <- (p2 + p4) +
  plot_annotation(title = "Pairwise Differences Between Drugs (Correct Trials)",
                  theme = theme(plot.title = element_text(face = "bold", size = fontsize)))

# Full combined plot
p_combined <- (p_top / p_bottom) +
  plot_layout(heights = c(3, 2))

# Save combined
combined_filename <- sprintf("metacognition_combined_%s.png", TIMESTAMP)
outfile_combined <- file.path(FIGS_DIR, combined_filename)
ggsave(outfile_combined, p_combined, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved combined plot to: ", outfile_combined)

# Save individual panels
ggsave(file.path(FIGS_DIR, sprintf("hrd_predicted_confidence_%s.png", TIMESTAMP)),
       p1, width = 6, height = 6, dpi = 300, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("hrd_pairwise_contrast_%s.png", TIMESTAMP)),
       p2, width = 6, height = 5, dpi = 300, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("rrst_predicted_confidence_%s.png", TIMESTAMP)),
       p3, width = 6, height = 6, dpi = 300, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("rrst_pairwise_contrast_%s.png", TIMESTAMP)),
       p4, width = 6, height = 5, dpi = 300, bg = "white")

# ============================================================
# HRD: Drug Main Effect + Interaction Detail Plot
# ============================================================
message("Generating HRD confidence detail plot...")

drug_cols <- c("placebo" = "#2E3FAE", "propranolol" = "#C13A8C", "bisoprolol" = "#1FA187")

# A) Marginal Drug means (averaging over ResponseCorrect)
emm_drug_marginal <- emmeans(hrd_model, ~ Drug, data = HRD_trl_data)
emm_drug_df <- as.data.frame(emm_drug_marginal)

# Back-transform from log-odds to probability scale for interpretability
emm_drug_df$pred <- plogis(emm_drug_df$emmean)
emm_drug_df$pred_lo <- plogis(emm_drug_df$asymp.LCL)
emm_drug_df$pred_hi <- plogis(emm_drug_df$asymp.UCL)

# Subject-level raw mean confidence (marginal over accuracy)
subj_means_marginal <- HRD_trl_data %>%
  group_by(subject, Drug) %>%
  summarise(mean_conf = mean(Conf, na.rm = TRUE), .groups = "drop")

p_hrd_marginal <- ggplot(emm_drug_df, aes(x = Drug, y = pred, color = Drug)) +
  geom_jitter(data = subj_means_marginal, aes(x = Drug, y = mean_conf),
              width = 0.15, alpha = 0.25, size = 1.5, show.legend = FALSE) +
  geom_pointrange(aes(ymin = pred_lo, ymax = pred_hi),
                  size = 1.2, linewidth = 1.2, show.legend = FALSE) +
  scale_color_manual(values = drug_cols) +
  labs(x = "", y = "Predicted Confidence",
       title = "Drug Main Effect") +
  coord_cartesian(ylim = c(0.3, 0.85)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold"),
    axis.line = element_line(linewidth = 0.8)
  )

# B) Drug x Accuracy cell means
emm_cell <- emmeans(hrd_model, ~ Drug | ResponseCorrect, data = HRD_trl_data)
emm_cell_df <- as.data.frame(emm_cell)
emm_cell_df$pred <- plogis(emm_cell_df$emmean)
emm_cell_df$pred_lo <- plogis(emm_cell_df$asymp.LCL)
emm_cell_df$pred_hi <- plogis(emm_cell_df$asymp.UCL)

# Subject-level raw means by Drug x Accuracy
subj_means_cell <- HRD_trl_data %>%
  group_by(subject, Drug, ResponseCorrect) %>%
  summarise(mean_conf = mean(Conf, na.rm = TRUE), .groups = "drop")

p_hrd_interaction <- ggplot(emm_cell_df, aes(x = Drug, y = pred, color = Drug)) +
  geom_jitter(data = subj_means_cell, aes(x = Drug, y = mean_conf),
              width = 0.15, alpha = 0.2, size = 1.2, show.legend = FALSE) +
  geom_pointrange(aes(ymin = pred_lo, ymax = pred_hi),
                  size = 1.0, linewidth = 1.0, show.legend = FALSE) +
  facet_wrap(~ ResponseCorrect) +
  scale_color_manual(values = drug_cols) +
  labs(x = "", y = "",
       title = "Drug \u00d7 Accuracy Interaction") +
  coord_cartesian(ylim = c(0.3, 0.85)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    axis.line = element_line(linewidth = 0.8)
  )

p_hrd_detail <- p_hrd_marginal + p_hrd_interaction +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title = "HRD Metacognitive Confidence: Drug Effects",
    subtitle = "Points = subject means; error bars = model-estimated marginal means \u00b1 95% CI",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "grey40")
    )
  )

ggsave(file.path(FIGS_DIR, sprintf("hrd_confidence_drug_detail_%s.png", TIMESTAMP)),
       p_hrd_detail, width = 14, height = 6, dpi = 300, bg = "white")
message("Saved HRD confidence detail plot.")

message("Plot generation complete!")
