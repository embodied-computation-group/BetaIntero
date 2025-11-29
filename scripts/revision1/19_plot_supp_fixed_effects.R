#!/usr/bin/env Rscript
# Supplementary plot: drug main effects, nTrials, and visit effects
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

# Load most recent model
model_files <- list.files(DATA_DIR, pattern = "^model_raw_interaction_.*\\.rds$", full.names = TRUE)
if (length(model_files) == 0) stop("No model files found in data directory: ", DATA_DIR)
model_file <- model_files[order(file.info(model_files)$mtime, decreasing = TRUE)][1]
message("Loading model from: ", model_file)
mod <- readRDS(model_file)

# Create plots
# 1) Drug main effect
p_drug <- plot_model(mod, type = "pred", terms = "drugs", title = "Drug: Main Effect", axis.title = c("Drug", "Pr(Response = More)")) + theme_sjplot() + labs(x = "Drug")

# 2) nTrials effect (continuous)
p_nTrials <- plot_model(mod, type = "pred", terms = "nTrials [all]", title = "Trial Number Effect", axis.title = c("Trial Number", "Pr(Response = More)")) + theme_sjplot() + labs(x = "Trial Number")

# 3) Visit effect (factor)
p_visit <- plot_model(mod, type = "pred", terms = "visit", title = "Visit Effect", axis.title = c("Visit", "Pr(Response = More)")) + theme_sjplot() + labs(x = "Visit")

# Combine
p_comb <- (p_drug | p_nTrials | p_visit) + plot_annotation(title = "Supplementary: Lesser Fixed Effects (Drug, Trial, Visit)")

# Save
png_file <- file.path(FIGS_DIR, "supp_fixed_effects.png")
svg_file <- file.path(FIGS_DIR, "supp_fixed_effects.svg")
ggsave(png_file, plot = p_comb, width = 14, height = 5, dpi = 300, bg = "white")
ggsave(svg_file, plot = p_comb, width = 14, height = 5)
message("Saved supplementary fixed-effects plots to:", png_file, " and ", svg_file)
