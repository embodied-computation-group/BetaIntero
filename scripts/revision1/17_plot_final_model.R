#!/usr/bin/env Rscript
# Plot main effects and drug * listenBPM interaction for the final model (no 3-way interaction)

suppressPackageStartupMessages({
  library(glmmTMB)
  library(ggplot2)
  library(sjPlot)
  library(here)
})

# Directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "revision1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)

# Load final model (no 3-way interaction) - choose the most recent saved model
model_files <- list.files(DATA_DIR, pattern = "^model_raw_interaction_.*\\.rds$", full.names = TRUE)
if (length(model_files) == 0) stop("No model files found in data directory: ", DATA_DIR)
model_file <- model_files[order(file.info(model_files)$mtime, decreasing = TRUE)][1]
message("Loading model from: ", model_file)
mod <- readRDS(model_file)

# Plot main effect: Listen BPM
p_listen_main <- plot_model(mod, type = "pred", terms = "listenBPM_centered [all]",
                           title = "Main Effect of Listen BPM",
                           axis.title = c("Listen BPM (Centered)", "Pr(Response = More)")) +
  theme_sjplot()

# Plot main effect: Response BPM
p_response_main <- plot_model(mod, type = "pred", terms = "responseBPM_centered [all]",
                            title = "Main Effect of Response BPM",
                            axis.title = c("Response BPM (Centered)", "Pr(Response = More)")) +
  theme_sjplot()

# Plot drug * listenBPM interaction
p_drug_listen <- plot_model(mod, type = "pred", terms = c("listenBPM_centered [all]", "drugs"),
                          title = "Interaction: Drug × Listen BPM",
                          axis.title = c("Listen BPM (Centered)", "Pr(Response = More)")) +
  theme_sjplot() + labs(color = "Drug")

# Combine plots
library(patchwork)
p_combined <- (p_listen_main | p_response_main) / p_drug_listen +
  plot_annotation(tag_levels = 'A', title = "Model Results: Raw BPM Analysis (Final Model)")

# Save combined plot
combined_filename <- "combined_raw_results_final_model.png"
combined_file <- file.path(FIGS_DIR, combined_filename)
ggsave(combined_file, plot = p_combined, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved combined plot to: ", combined_file)

# Save individual plots
ggsave(file.path(FIGS_DIR, "main_effect_listenBPM_final.png"), p_listen_main, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, "main_effect_responseBPM_final.png"), p_response_main, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, "interaction_drug_listenBPM_final.png"), p_drug_listen, width = 8, height = 6, bg = "white")
message("Plotting complete!")
