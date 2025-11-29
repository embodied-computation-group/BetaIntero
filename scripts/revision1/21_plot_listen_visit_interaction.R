#!/usr/bin/env Rscript
# Plot Listen BPM × Visit interaction from the listen-visit model
suppressPackageStartupMessages({
  library(glmmTMB)
  library(ggplot2)
  library(sjPlot)
  library(patchwork)
  library(here)
})

# Directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "revision1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)

# Find most recent listen-visit model
model_files <- list.files(DATA_DIR, pattern = "^model_listen_visit_.*\\.rds$", full.names = TRUE)
if (length(model_files) == 0) stop("No listen-visit model files found in data directory: ", DATA_DIR)
model_file <- model_files[order(file.info(model_files)$mtime, decreasing = TRUE)][1]
message("Loading model from: ", model_file)
mod <- readRDS(model_file)

# Plot the interaction using sjPlot (bias_correction recommended for GLMMs)
# Show predicted probability curves for Listen BPM across Visit levels
p_inter <- tryCatch({
  plot_model(mod, type = "pred", terms = c("listenBPM_centered [all]", "visit"),
             title = "Interaction: Listen BPM × Visit",
             axis.title = c("Listen BPM (Centered)", "Pr(Response = More)"),
             legend.title = "Visit",
             show.data = FALSE,
             bias_correction = TRUE)
}, error = function(e) {
  message("plot_model failed with bias_correction=TRUE; retrying without bias correction: ", e$message)
  plot_model(mod, type = "pred", terms = c("listenBPM_centered [all]", "visit"),
             title = "Interaction: Listen BPM × Visit",
             axis.title = c("Listen BPM (Centered)", "Pr(Response = More)"),
             legend.title = "Visit",
             show.data = FALSE)
})

# Tidy up theme
p_inter <- p_inter + theme_minimal() + theme(legend.position = "bottom")

# Save
png_file <- file.path(FIGS_DIR, "listen_visit_interaction.png")
svg_file <- file.path(FIGS_DIR, "listen_visit_interaction.svg")
ggsave(png_file, plot = p_inter, width = 8, height = 5, dpi = 300, bg = "white")
ggsave(svg_file, plot = p_inter, width = 8, height = 5)
message("Saved Listen×Visit interaction plots to:\n", png_file, "\n", svg_file)
