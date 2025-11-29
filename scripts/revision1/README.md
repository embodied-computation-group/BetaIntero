# Revision 1 Analysis

This directory contains scripts for the revision 1 analysis of the BetaIntero project.

## Scripts

- `01_intero_interaction_model.R`: Fits hierarchical binomial logistic regression models predicting choice (Decision) from `listenBPM`, `responseBPM`, and their interaction.
  - Models are fitted separately for each drug condition (PLACEBO, BISO, PROP).
  - Uses `glmmTMB` with maximal random slopes structure where possible, degrading gracefully to simpler structures if convergence fails.
  - Outputs are saved in `outputs/`.

## Outputs

The `outputs/` directory (created by the script) will contain:
- `figs/`: Figures (currently empty, script focuses on models).
- `data/`: Saved model objects (`.rds`) and summary tables (`.csv`).
- `docs/`: Text summaries of the models.
