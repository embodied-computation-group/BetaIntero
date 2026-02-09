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
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
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
emm_hrd <- emmeans(hrd_model, ~ Drug | ResponseCorrect)
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
emm_rrst <- emmeans(rrst_model, ~ Drug | Accuracy)
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

message("Plot generation complete!")
