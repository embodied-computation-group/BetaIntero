#!/usr/bin/env Rscript
# Simple 3-way interaction plot for Drug × Listen BPM × Response BPM
# Panelled by drug, lines for listenBPM_centered and responseBPM_centered

if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org")
pacman::p_load(glmmTMB, sjPlot, here)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "revision1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")

dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)

# Load most recent model
model_pattern <- "model_raw_interaction_.*\\.rds"
model_files <- list.files(DATA_DIR, pattern = model_pattern, full.names = TRUE)
if (length(model_files) == 0) stop("No raw model files found in ", DATA_DIR)
details <- file.info(model_files)
model_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
mod <- readRDS(model_file)

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Plot: lines for listenBPM_centered and responseBPM_centered, panelled by drug

# Panel 1: lines for listenBPM_centered, panelled by drug
p_listen <- plot_model(
  mod,
  type = "pred",
  terms = c("drugs", "listenBPM_centered [all]"),
  title = "Effect of Listen BPM (Centered) by Drug",
  axis.title = c("Listen BPM (Centered)", "Pr(Response = More)")
) + theme_sjplot()

# Panel 2: lines for responseBPM_centered, panelled by drug
p_response <- plot_model(
  mod,
  type = "pred",
  terms = c("drugs", "responseBPM_centered [all]"),
  title = "Effect of Response BPM (Centered) by Drug",
  axis.title = c("Response BPM (Centered)", "Pr(Response = More)")
) + theme_sjplot()

# Save plots
listen_file <- file.path(FIGS_DIR, sprintf("listenBPM_by_drug_%s.png", TIMESTAMP))
response_file <- file.path(FIGS_DIR, sprintf("responseBPM_by_drug_%s.png", TIMESTAMP))
ggsave(listen_file, plot = p_listen, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(response_file, plot = p_response, width = 8, height = 6, dpi = 300, bg = "white")

message("Saved listenBPM plot to: ", listen_file)
message("Saved responseBPM plot to: ", response_file)
message("Plotting complete!")
