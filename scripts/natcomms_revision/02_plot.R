#!/usr/bin/env Rscript
# Plot results for the combined drug interaction model with Trial Number control

# Ensure required packages are installed
if (!require("pacman")) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
pacman::p_load(dplyr, readr, glmmTMB, ggplot2, sjPlot, ggeffects, here, patchwork)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")

# Load model
# Find most recent model file matching the pattern
model_pattern <- "model_combined_drug_interaction_with_trialN_.*\\.rds"
model_files <- list.files(DATA_DIR, pattern = model_pattern, full.names = TRUE)

if (length(model_files) == 0) {
  # Fallback to legacy name if no timestamped files found
  legacy_file <- file.path(DATA_DIR, "model_combined_drug_interaction_with_trialN.rds")
  if (file.exists(legacy_file)) {
    model_file <- legacy_file
  } else {
    stop("No model files found in ", DATA_DIR)
  }
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

# 1. Fixed Effects Forest Plot
# User requested highlighting the OR=1 line and clipping the x-axis at 100
p_forest <- plot_model(mod, type = "est", 
                       show.values = TRUE, show.p = TRUE,
                       title = "Fixed Effects",
                       vline.color = "red",   # Highlight the null line (OR = 1)
                       ylim = c(0.1, 100)) +  # Clip x-axis (OR) to [0.1, 100]
  theme_sjplot()

# 2. Marginal Effects Plots

# A. Listen BPM (Main Effect - averaged over drugs)
p_listen <- plot_model(mod, type = "pred", terms = "listenBPM_centered [all]",
                       title = "Effect of Listen BPM",
                       axis.title = c("Listen BPM (BPM, centered)", "Pr(Response = More)")) +
  theme_sjplot()

# B. Response BPM (Main Effect)
p_response <- plot_model(mod, type = "pred", terms = "responseBPM_centered [all]",
                         title = "Effect of Response BPM",
                         axis.title = c("Response BPM (BPM, centered)", "Pr(Response = More)")) +
  theme_sjplot()

# C. Trial Number (Control Effect)
p_trial <- plot_model(mod, type = "pred", terms = "nTrials [all]",
                   title = "Effect of Trial Number",
                   axis.title = c("Trial Number", "Pr(Response = More)")) +
  theme_sjplot()

# D. Visit (Control Effect)
p_visit <- plot_model(mod, type = "pred", terms = "visit",
                      title = "Effect of Visit",
                      axis.title = c("Visit", "Pr(Response = More)")) +
  theme_sjplot()

# E. Drug * Listen BPM Interaction
p_int <- plot_model(mod, type = "pred", terms = c("listenBPM_centered [all]", "drugs"),
                    title = "Interaction: Drug * Listen BPM",
                    axis.title = c("Listen BPM (BPM, centered)", "Pr(Response = More)")) +
  theme_sjplot() +
  labs(color = "Drug")

# Combine plots using patchwork
# Layout: 
# Row 1: Forest Plot
# Row 2: Listen BPM | Response BPM | Trial Number
# Row 3: Visit | Interaction
p_combined <- (p_forest) / (p_listen + p_response + p_trial) / (p_visit + p_int) +
  plot_annotation(tag_levels = 'A', title = "Model Results: Interoception Choice (with Trial N & Visit Control)")

# Save combined plot
combined_filename <- sprintf("combined_model_results_with_trialN_%s.png", TIMESTAMP)
outfile_combined <- file.path(FIGS_DIR, combined_filename)
ggsave(outfile_combined, p_combined, width = 12, height = 12, bg = "white")
message("Saved combined plot to: ", outfile_combined)

# Save individual plots as well
ggsave(file.path(FIGS_DIR, sprintf("forest_plot_trialN_%s.png", TIMESTAMP)), p_forest, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("effect_listenBPM_trialN_%s.png", TIMESTAMP)), p_listen, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("effect_responseBPM_trialN_%s.png", TIMESTAMP)), p_response, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("effect_trialN_%s.png", TIMESTAMP)), p_trial, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("effect_visit_%s.png", TIMESTAMP)), p_visit, width = 6, height = 5, bg = "white")
ggsave(file.path(FIGS_DIR, sprintf("interaction_drug_listenBPM_trialN_%s.png", TIMESTAMP)), p_int, width = 7, height = 5, bg = "white")

# 3. Combo plot: A) Response BPM, B) Listen BPM, C) Drug * Listen BPM Interaction
# Layout: (A + B) / C
p_combo <- (p_response + p_listen) / p_int +
  plot_annotation(tag_levels = 'A')

combo_filename <- sprintf("combo_main_effects_%s.png", TIMESTAMP)
outfile_combo <- file.path(FIGS_DIR, combo_filename)
ggsave(outfile_combo, p_combo, width = 10, height = 8, bg = "white")
message("Saved combo plot to: ", outfile_combo)

# --- 4. Journal-style combo plot ---
# Colour palette
drug_cols <- c("PLACEBO" = "#2E3FAE", "BISO" = "#1FA187", "PROP" = "#C13A8C")
drug_labels <- c("PLACEBO" = "Placebo", "BISO" = "Bisoprolol", "PROP" = "Propranolol")

# Journal theme: white background, no gridlines, black axes/text
theme_journal <- theme_classic(base_size = 12) +
  theme(
    text             = element_text(colour = "black"),
    axis.text        = element_text(colour = "black"),
    axis.ticks       = element_line(colour = "black"),
    axis.line        = element_line(colour = "black"),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0),
    legend.title     = element_text(face = "bold"),
    legend.position  = "bottom",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    strip.background = element_blank()
  )

# Use ggpredict so we have full control over the ggplot layers
pred_response <- ggpredict(mod, terms = "responseBPM_centered [all]")
pred_listen   <- ggpredict(mod, terms = "listenBPM_centered [all]")
pred_int      <- ggpredict(mod, terms = c("listenBPM_centered [all]", "drugs"))

# A) Response BPM main effect
pj_response <- ggplot(pred_response, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.8, colour = "black") +
  labs(x = "Response BPM (centered)", y = "Pr(Response = More)") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_journal

# B) Listen BPM main effect
pj_listen <- ggplot(pred_listen, aes(x, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.8, colour = "black") +
  labs(x = "Listen BPM (centered)", y = "Pr(Response = More)") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_journal

# C) Drug * Listen BPM interaction
pj_int <- ggplot(pred_int, aes(x, predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = drug_cols, labels = drug_labels) +
  scale_fill_manual(values = drug_cols, labels = drug_labels) +
  labs(x = "Listen BPM (centered)", y = "Pr(Response = More)",
       colour = "Drug", fill = "Drug") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_journal

# Assemble: (A + B) / C
pj_combo <- (pj_response + pj_listen) / pj_int +
  plot_annotation(tag_levels = 'A',
                  theme = theme(plot.tag = element_text(face = "bold", size = 14)))

combo_j_filename <- sprintf("combo_journal_%s.png", TIMESTAMP)
outfile_combo_j <- file.path(FIGS_DIR, combo_j_filename)
ggsave(outfile_combo_j, pj_combo, width = 10, height = 8, dpi = 300, bg = "white")
message("Saved journal combo plot to: ", outfile_combo_j)

message("Analysis complete!")
