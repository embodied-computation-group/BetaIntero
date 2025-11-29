#!/usr/bin/env Rscript
# Plot 2-way interaction: listenBPM_centered by responseBPM_centered
# No drug panels, just the interaction surface

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

# Plot: 2-way interaction listenBPM_centered by responseBPM_centered
p_2way <- plot_model(
  mod,
  type = "pred",
  terms = c("listenBPM_centered [all]", "responseBPM_centered [all]"),
  title = "Interaction: Listen BPM × Response BPM",
  axis.title = c("Listen BPM (Centered)", "Pr(Response = More)")
) + theme_sjplot()

# Save plot
plot_file <- file.path(FIGS_DIR, sprintf("listen_by_response_interaction_%s.png", TIMESTAMP))
ggsave(plot_file, plot = p_2way, width = 8, height = 6, dpi = 300, bg = "white")

message("Saved 2-way interaction plot to: ", plot_file)
message("Plotting complete!")
