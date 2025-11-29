#!/usr/bin/env Rscript
# Plot main effect of drug on Pr(Response = More)

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

# Plot main effect of drug
p_drug <- plot_model(
  mod,
  type = "pred",
  terms = "drugs",
  title = "Main Effect of Drug on Pr(Response = More)",
  axis.title = c("Drug", "Pr(Response = More)")
) + theme_sjplot()

outfile_drug <- file.path(FIGS_DIR, "main_effect_drug.png")
ggsave(outfile_drug, p_drug, width = 6, height = 5, bg = "white")
message("Saved main effect plot to: ", outfile_drug)
