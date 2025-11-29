#!/usr/bin/env Rscript
# Plot results for the Raw BPM model
# Model: Response ~ drugs * listenBPM + responseBPM + nTrials + visit

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
model_pattern <- "model_raw_interaction_.*\\.rds"
model_files <- list.files(DATA_DIR, pattern = model_pattern, full.names = TRUE)

if (length(model_files) == 0) {
  stop("No raw model files found in ", DATA_DIR)
} else {
  # Sort by modification time to get the most recent
  details <- file.info(model_files)
  model_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
}

message("Loading model from: ", basename(model_file))
mod <- readRDS(model_file)

# Generate Timestamp for plots
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Plot Run Timestamp: ", TIMESTAMP)

message("Generating plots...")


# --- 3-way Interaction Plot: Drug × Listen BPM × Response BPM ---
# Use sjPlot to visualize the 3-way interaction as predicted probabilities
p_3way <- plot_model(
  mod,
  type = "pred",
  terms = c("listenBPM_centered [all]", "responseBPM_centered [all]", "drugs"),
  title = "3-way Interaction: Drug × Listen BPM × Response BPM",
  axis.title = c("Listen BPM (Centered)", "Pr(Response = More)")
) + theme_sjplot() + labs(color = "Drug")
# --- 1. Forest Plot of Fixed Effects ---
p_forest <- plot_model(mod, type = "est", 
                       show.values = TRUE, show.p = TRUE,
                       title = "Fixed Effects (Raw BPM Model)",
                       vline.color = "red",
                       ylim = c(0.1, 100)) + # Clip large ORs if necessary
  theme_sjplot()

# --- 2. Interaction Plot: Drug * Listen BPM ---
# We want to see how the probability of "More" changes with Listen BPM for each drug
# holding responseBPM constant (at mean).
p_listen <- plot_model(mod, type = "pred", terms = c("listenBPM_centered [all]", "drugs"),
                       title = "Interaction: Drug * Listen BPM",
                       axis.title = c("Listen BPM (Centered)", "Pr(Response = More)")) +
  theme_sjplot() +
  labs(color = "Drug")

# --- 3. Main Effect: Response BPM ---
p_response <- plot_model(mod, type = "pred", terms = "responseBPM_centered [all]",
                         title = "Effect of Response BPM",
                         axis.title = c("Response BPM (Centered)", "Pr(Response = More)")) +
  theme_sjplot()


# Save individual 3-way plot first to ensure it is always written
ggsave(file.path(FIGS_DIR, sprintf("interaction_3way_raw_%s.png", TIMESTAMP)), p_3way, width = 10, height = 5, bg = "white")

# --- 5. Combine Plots ---
p_combined <- (p_forest | p_response) / p_listen / p_3way +
  plot_annotation(tag_levels = 'A', title = "Model Results: Raw BPM Analysis (3-way Interaction)")

# Save combined plot
combined_filename <- sprintf("combined_raw_results_3way_%s.png", TIMESTAMP)
combined_file <- file.path(FIGS_DIR, combined_filename)
ggsave(combined_file, plot = p_combined, width = 14, height = 14, dpi = 300, bg = "white")
message("Saved combined plot to: ", combined_file)

# Save individual plots
ggsave(file.path(FIGS_DIR, sprintf("forest_plot_raw_%s.png", TIMESTAMP)), p_forest, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("interaction_listenBPM_raw_%s.png", TIMESTAMP)), p_listen, width = 8, height = 6, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("effect_responseBPM_raw_%s.png", TIMESTAMP)), p_response, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("interaction_3way_raw_%s.png", TIMESTAMP)), p_3way, width = 10, height = 5, bg = "white")

message("Plotting complete!")
