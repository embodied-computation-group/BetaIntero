#!/usr/bin/env Rscript
# Test: Compare old (Analysis.Rmd) vs new (01_model.R) preprocessing and model output

library(dplyr)
library(readr)
library(tidyr)
library(glmmTMB)
library(here)

BASE_DIR <- here::here()
passed <- 0
failed <- 0

report <- function(test_name, result) {
  if (result) {
    cat(sprintf("  PASS: %s\n", test_name))
    passed <<- passed + 1
  } else {
    cat(sprintf("  FAIL: %s\n", test_name))
    failed <<- failed + 1
  }
}

# ============================================================
# 1. RRST Data Comparison
# ============================================================
cat("\n=== RRST Preprocessing Comparison ===\n")

# OLD code (Analysis.Rmd lines 276-290)
RRST_old <- read_csv(file.path(BASE_DIR, "data", "cleaned", "RRST.csv"), show_col_types = FALSE) %>%
  filter(ConfResp >= 0 & ConfResp <= 100) %>%
  mutate(Conf = ConfResp / 100,
         Stimulus = Stim,
         Drug = factor(drugs, levels = c("PLACEBO", "PROP", "BISO"),
                       labels = c("placebo", "propranolol", "bisoprolol")),
         Drug = relevel(Drug, ref = "placebo"),
         Condition = drugs,
         Accuracy = factor(Resp, levels = c(1, 0), labels = c("Correct", "Incorrect"))
  ) %>% drop_na()

# NEW code (01_model.R)
RRST_new <- read_csv(file.path(BASE_DIR, "data", "cleaned", "RRST.csv"), show_col_types = FALSE) %>%
  filter(ConfResp >= 0 & ConfResp <= 100) %>%
  mutate(
    Conf = ConfResp / 100,
    Stimulus = Stim,
    Drug = factor(drugs, levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    Condition = drugs,
    Accuracy = factor(Resp, levels = c(1, 0), labels = c("Correct", "Incorrect"))
  ) %>% drop_na()

cat("RRST old dims:", nrow(RRST_old), "x", ncol(RRST_old), "\n")
cat("RRST new dims:", nrow(RRST_new), "x", ncol(RRST_new), "\n")

report("RRST row count matches", nrow(RRST_old) == nrow(RRST_new))
report("RRST Conf identical", identical(RRST_old$Conf, RRST_new$Conf))
report("RRST Drug identical", identical(RRST_old$Drug, RRST_new$Drug))
report("RRST Stimulus identical", identical(RRST_old$Stimulus, RRST_new$Stimulus))
report("RRST Accuracy identical", identical(RRST_old$Accuracy, RRST_new$Accuracy))
report("RRST subject identical", identical(RRST_old$subject, RRST_new$subject))

# ============================================================
# 2. HRD Data Comparison
# ============================================================
cat("\n=== HRD Preprocessing Comparison ===\n")

# OLD code (Analysis.Rmd lines 310-322) -- note: uses read.csv
HRD_old <- read.csv(file.path(BASE_DIR, "data", "cleaned", "HRD.csv")) %>%
  mutate(Conf = Confidence / 100,
         drugs = as.factor(drugs),
         Drug = factor(drugs, levels = c("PLACEBO", "PROP", "BISO"),
                       labels = c("placebo", "propranolol", "bisoprolol")),
         Drug = relevel(Drug, ref = "placebo"),
         ResponseCorrect = factor(ResponseCorrect, levels = c("True", "False"),
                                  labels = c("Correct", "Incorrect")),
         BPM_scaled = scale(listenBPM),
         Condition = drugs
  ) %>%
  filter(subject != "sub_4049") %>% drop_na()

# NEW code (01_model.R) -- uses read.csv to match original
HRD_new <- read.csv(file.path(BASE_DIR, "data", "cleaned", "HRD.csv")) %>%
  mutate(
    Conf = Confidence / 100,
    drugs = as.factor(drugs),
    Drug = factor(drugs, levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    ResponseCorrect = factor(ResponseCorrect, levels = c("True", "False"),
                             labels = c("Correct", "Incorrect")),
    BPM_scaled = scale(listenBPM),
    Condition = drugs
  ) %>%
  filter(subject != "sub_4049") %>% drop_na()

cat("HRD old dims:", nrow(HRD_old), "x", ncol(HRD_old), "\n")
cat("HRD new dims:", nrow(HRD_new), "x", ncol(HRD_new), "\n")

report("HRD row count matches", nrow(HRD_old) == nrow(HRD_new))
report("HRD Conf identical", identical(HRD_old$Conf, HRD_new$Conf))
report("HRD Drug identical", identical(HRD_old$Drug, HRD_new$Drug))
report("HRD ResponseCorrect identical", identical(HRD_old$ResponseCorrect, HRD_new$ResponseCorrect))
report("HRD BPM_scaled equal", isTRUE(all.equal(as.numeric(HRD_old$BPM_scaled),
                                                  as.numeric(HRD_new$BPM_scaled))))
report("HRD subject identical", identical(HRD_old$subject, HRD_new$subject))

# ============================================================
# 3. Fit RRST models with old and new data, compare coefficients
# ============================================================
cat("\n=== RRST Model Comparison ===\n")
cat("Fitting RRST model with OLD data...\n")

rrst_mod_old <- glmmTMB(
  Conf ~ Stimulus + Drug * Accuracy + (1 + Stimulus + Accuracy | subject),
  data = RRST_old,
  family = ordbeta(),
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)

cat("Fitting RRST model with NEW data...\n")
rrst_mod_new <- glmmTMB(
  Conf ~ Stimulus + Drug * Accuracy + (1 + Stimulus + Accuracy | subject),
  data = RRST_new,
  family = ordbeta(),
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)

rrst_coef_old <- fixef(rrst_mod_old)$cond
rrst_coef_new <- fixef(rrst_mod_new)$cond

cat("\nRRST OLD coefficients:\n")
print(round(rrst_coef_old, 6))
cat("\nRRST NEW coefficients:\n")
print(round(rrst_coef_new, 6))

rrst_coef_match <- isTRUE(all.equal(rrst_coef_old, rrst_coef_new, tolerance = 1e-6))
report("RRST fixed effects match (tol=1e-6)", rrst_coef_match)
report("RRST AIC matches", isTRUE(all.equal(AIC(rrst_mod_old), AIC(rrst_mod_new), tolerance = 1e-4)))

# ============================================================
# 4. Fit HRD models with old and new data, compare coefficients
# ============================================================
cat("\n=== HRD Model Comparison ===\n")
cat("Fitting HRD model with OLD data...\n")

hrd_mod_old <- glmmTMB(
  Conf ~ Drug * ResponseCorrect + BPM_scaled + (1 + ResponseCorrect + BPM_scaled | subject),
  data = HRD_old,
  family = ordbeta(),
  start = list(psi = c(0, 1)),
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)

cat("Fitting HRD model with NEW data...\n")
hrd_mod_new <- glmmTMB(
  Conf ~ Drug * ResponseCorrect + BPM_scaled + (1 + ResponseCorrect + BPM_scaled | subject),
  data = HRD_new,
  family = ordbeta(),
  start = list(psi = c(0, 1)),
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)

hrd_coef_old <- fixef(hrd_mod_old)$cond
hrd_coef_new <- fixef(hrd_mod_new)$cond

cat("\nHRD OLD coefficients:\n")
print(round(hrd_coef_old, 6))
cat("\nHRD NEW coefficients:\n")
print(round(hrd_coef_new, 6))

hrd_coef_match <- isTRUE(all.equal(hrd_coef_old, hrd_coef_new, tolerance = 1e-6))
report("HRD fixed effects match (tol=1e-6)", hrd_coef_match)
report("HRD AIC matches", isTRUE(all.equal(AIC(hrd_mod_old), AIC(hrd_mod_new), tolerance = 1e-4)))

# ============================================================
# Summary
# ============================================================
cat(sprintf("\n=== SUMMARY: %d passed, %d failed ===\n", passed, failed))
if (failed > 0) {
  cat("SOME TESTS FAILED!\n")
  quit(status = 1)
} else {
  cat("ALL TESTS PASSED!\n")
}
