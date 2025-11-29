#!/usr/bin/env Rscript
# Plot listenBPM_centered × responseBPM_centered interaction for latest model

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

# Load most recent model
model_pattern <- "model_raw_interaction_.*\\.rds"
model_files <- list.files(DATA_DIR, pattern = model_pattern, full.names = TRUE)
stopifnot(length(model_files) > 0)
details <- file.info(model_files)
model_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
mod <- readRDS(model_file)

# Plot interaction: x = responseBPM_centered, lines = listenBPM_centered at -20, 0, +20
p_int <- plot_model(
  mod,
  type = "pred",
  terms = c("responseBPM_centered [all]", "listenBPM_centered [-20,0,20]"),
  title = "Interaction: Response BPM × Listen BPM (Centered)",
  axis.title = c("Response BPM (Centered)", "Pr(Response = More)")
) + theme_sjplot() + labs(color = "Listen BPM\n(Centered)")

outfile_int <- file.path(FIGS_DIR, "interaction_plot_listen_by_response.png")
ggsave(outfile_int, p_int, width = 8, height = 6, bg = "white")
message("Saved plot to: ", outfile_int)
