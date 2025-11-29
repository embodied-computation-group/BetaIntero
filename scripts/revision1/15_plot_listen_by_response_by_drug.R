#!/usr/bin/env Rscript
# Plot listenBPM_centered × responseBPM_centered interaction, panelled by drug

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




# Plot interaction: x = listenBPM_centered, lines = responseBPM_centered at -15, -10, -5, 0, 5, 10, 15, panelled by drug
p_int_panel <- plot_model(
  mod,
  type = "pred",
  terms = c("listenBPM_centered [all]", "responseBPM_centered [-15,-10,-5,0,5,10,15]", "drugs"),
  title = "Interaction: Listen BPM × Response BPM (Centered) by Drug",
  axis.title = c("Listen BPM (Centered)", "Pr(Response = More)")
) + theme_sjplot() + labs(color = "Response BPM\n(Centered)")

outfile_panel <- file.path(FIGS_DIR, "interaction_plot_listen_by_response_by_drug.png")
ggsave(outfile_panel, p_int_panel, width = 12, height = 6, bg = "white")
message("Saved panelled plot to: ", outfile_panel)
