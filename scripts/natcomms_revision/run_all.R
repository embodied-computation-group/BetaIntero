#!/usr/bin/env Rscript
# Run all Nature Communications revision analyses
# 1. HRD logit (interoceptive choice)
# 2. Metacognition core (RRST + HRD confidence)
# 3. Metacognition RRST interaction (Stimulus x Drug x Accuracy)

library(here)

BASE_DIR <- here::here()
LOGIT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision", "hrd_logit")
META_DIR  <- file.path(BASE_DIR, "scripts", "natcomms_revision", "metacognition")

run_script <- function(script_path) {
  label <- basename(dirname(script_path))
  name  <- basename(script_path)
  message("\n", strrep("=", 60))
  message(sprintf("  [%s] %s", label, name))
  message(strrep("=", 60), "\n")
  source(script_path, local = new.env(parent = globalenv()))
  message(sprintf("\n  Done: %s/%s\n", label, name))
}

t0 <- Sys.time()

# --- HRD Logit ---
run_script(file.path(LOGIT_DIR, "01_model.R"))
run_script(file.path(LOGIT_DIR, "02_plot.R"))
run_script(file.path(LOGIT_DIR, "03_diagnostics.R"))
run_script(file.path(LOGIT_DIR, "04_apa_table.R"))

# --- Metacognition (core RRST + HRD) ---
run_script(file.path(META_DIR, "01_model.R"))
run_script(file.path(META_DIR, "02_plot.R"))
run_script(file.path(META_DIR, "03_diagnostics.R"))
run_script(file.path(META_DIR, "04_apa_table.R"))

# --- Metacognition (RRST interaction) ---
run_script(file.path(META_DIR, "05_rrst_interaction.R"))
run_script(file.path(META_DIR, "06_rrst_interaction_diagnostics.R"))
run_script(file.path(META_DIR, "07_rrst_interaction_apa_table.R"))

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
message("\n", strrep("=", 60))
message("  All analyses complete. Total time: ", elapsed, " min")
message(strrep("=", 60))
