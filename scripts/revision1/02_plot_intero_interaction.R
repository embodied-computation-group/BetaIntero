#!/usr/bin/env Rscript
# Plot results for the combined drug interaction model

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

# Load model
model_file <- file.path(DATA_DIR, "model_combined_drug_interaction_with_HR.rds")
if (!file.exists(model_file)) stop("Model file not found: ", model_file)

message("Loading model from: ", model_file)
mod <- readRDS(model_file)

message("Generating plots...")

# 1. Fixed Effects Forest Plot
p_forest <- plot_model(mod, type = "est", 
                       show.values = TRUE, show.p = TRUE,
                       title = "Fixed Effects") +
  theme_sjplot()

# 2. Marginal Effects Plots

# A. Listen BPM (Main Effect - averaged over drugs)
p_listen <- plot_model(mod, type = "pred", terms = "listenBPM_scaled [all]",
                       title = "Effect of Listen BPM",
                       axis.title = c("Listen BPM (Scaled)", "Pr(Response = More)")) +
  theme_sjplot()

# B. Response BPM (Main Effect)
p_response <- plot_model(mod, type = "pred", terms = "responseBPM_scaled [all]",
                         title = "Effect of Response BPM",
                         axis.title = c("Response BPM (Scaled)", "Pr(Response = More)")) +
  theme_sjplot()

# C. Mean HR (Control Effect)
p_hr <- plot_model(mod, type = "pred", terms = "meanHR_scaled [all]",
                   title = "Effect of Mean HR",
                   axis.title = c("Mean HR (Scaled)", "Pr(Response = More)")) +
  theme_sjplot()

# D. Drug * Listen BPM Interaction
p_int <- plot_model(mod, type = "pred", terms = c("listenBPM_scaled [all]", "drugs"),
                    title = "Interaction: Drug * Listen BPM",
                    axis.title = c("Listen BPM (Scaled)", "Pr(Response = More)")) +
  theme_sjplot() +
  labs(color = "Drug")

# Combine plots using patchwork
# Layout: 
# Row 1: Forest Plot
# Row 2: Listen BPM | Response BPM | Mean HR
# Row 3: Interaction
p_combined <- (p_forest) / (p_listen + p_response + p_hr) / (p_int) +
  plot_annotation(tag_levels = 'A', title = "Model Results: Interoception Choice (with HR Control)")

# Save combined plot
outfile_combined <- file.path(FIGS_DIR, "combined_model_results_with_HR.png")
ggsave(outfile_combined, p_combined, width = 12, height = 12, bg = "white")
message("Saved combined plot to: ", outfile_combined)

# Save individual plots as well
ggsave(file.path(FIGS_DIR, "forest_plot.png"), p_forest, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, "effect_listenBPM.png"), p_listen, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, "effect_responseBPM.png"), p_response, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, "effect_meanHR.png"), p_hr, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, "interaction_drug_listenBPM.png"), p_int, width = 7, height = 5, bg = "white")

message("Analysis complete!")
