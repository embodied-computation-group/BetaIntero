#!/usr/bin/env Rscript
# Plot results for the Delta BPM model
# Model: Response ~ drugs * deltaBPM

# Ensure required packages are installed
if (!require("pacman")) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(dplyr, readr, glmmTMB, ggplot2, sjPlot, ggeffects, here, patchwork)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "revision1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")

# Create figs directory if it doesn't exist
dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)

# Load model
# Find most recent model file matching the pattern
model_pattern <- "model_diffBPM_interaction_.*\\.rds"
model_files <- list.files(DATA_DIR, pattern = model_pattern, full.names = TRUE)

if (length(model_files) == 0) {
  stop("No diffBPM model files found in ", DATA_DIR)
} else {
  # Sort by modification time to get the most recent
  details <- file.info(model_files)
  model_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
}

message("Loading model from: ", basename(model_file))
mod <- readRDS(model_file)

# Load PSE results
pse_pattern <- "pse_thresholds_diffBPM_.*\\.csv"
pse_files <- list.files(DATA_DIR, pattern = pse_pattern, full.names = TRUE)
if (length(pse_files) > 0) {
  details <- file.info(pse_files)
  pse_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
  message("Loading PSE results from: ", basename(pse_file))
  pse_df <- readr::read_csv(pse_file, show_col_types = FALSE)
} else {
  warning("No PSE results found.")
  pse_df <- NULL
}

# Generate Timestamp for plots
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Plot Run Timestamp: ", TIMESTAMP)

message("Generating plots...")

# --- 1. Forest Plot of Fixed Effects ---
p_forest <- plot_model(mod, type = "est", 
                       show.values = TRUE, show.p = TRUE,
                       title = "Fixed Effects (Diff BPM Model)",
                       vline.color = "red") +
  theme_sjplot()

# --- 2. Psychometric Curves (Predicted Probabilities) ---
# Effect of Diff BPM by Drug
p_curves <- plot_model(mod, type = "pred", terms = c("diffBPM [all]", "drugs"),
                    title = "Psychometric Curves: Drug * Diff BPM",
                    axis.title = c("Diff BPM (Response - Listen)", "Pr(Response = More)")) +
  theme_sjplot() +
  labs(color = "Drug") +
  geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5)

# --- 3. Threshold (PSE) Plot ---
if (!is.null(pse_df)) {
  # Ensure factor order
  pse_df$drugs <- factor(pse_df$drugs, levels = c("PLACEBO", "BISO", "PROP"))
  
  p_pse <- ggplot(pse_df, aes(x = drugs, y = PSE, fill = drugs)) +
    geom_bar(stat = "identity", alpha = 0.7, width = 0.6) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_sjplot() +
    labs(title = "Point of Subjective Equality (PSE)",
         y = "Threshold (Diff BPM at p=0.5)",
         x = "Drug Condition") +
    theme(legend.position = "none")
} else {
  p_pse <- ggplot() + theme_void() + labs(title = "PSE Data Missing")
}

# --- 4. Combine Plots ---
p_combined <- (p_forest | p_pse) / p_curves +
  plot_annotation(tag_levels = 'A', title = "Model Results: Diff BPM Analysis")

# Save combined plot
combined_filename <- sprintf("combined_diffBPM_results_%s.png", TIMESTAMP)
combined_file <- file.path(FIGS_DIR, combined_filename)
ggsave(combined_file, plot = p_combined, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved combined plot to: ", combined_file)

# Save individual plots
ggsave(file.path(FIGS_DIR, sprintf("forest_plot_diffBPM_%s.png", TIMESTAMP)), p_forest, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("psychometric_curves_diffBPM_%s.png", TIMESTAMP)), p_curves, width = 8, height = 6, bg = "white")
if (!is.null(pse_df)) {
  ggsave(file.path(FIGS_DIR, sprintf("pse_plot_diffBPM_%s.png", TIMESTAMP)), p_pse, width = 6, height = 5, bg = "white")
}

message("Plotting complete!")
