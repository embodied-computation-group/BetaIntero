#!/usr/bin/env Rscript
# Diagnostics for Metacognitive Confidence Models (RRST + HRD)
# 1. Confidence distribution (boundary check for ordbeta)
# 2. DHARMa residual diagnostics
# 3. VIF (collinearity) check
# 4. Random effects inspection and outlier detection

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, ggplot2, glmmTMB, here, tidyr, patchwork, performance, DHARMa)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision", "metacognition")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
FIGS_DIR <- file.path(OUTPUT_DIR, "figs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Diagnostics Run Timestamp: ", TIMESTAMP)

# --- Helper: clean_and_report ---
clean_and_report <- function(data, dataset_name) {
  total_rows_before <- nrow(data)
  cleaned_data <- data %>% tidyr::drop_na()
  total_rows_after <- nrow(cleaned_data)
  dropped_pct <- ((total_rows_before - total_rows_after) / total_rows_before) * 100
  message(sprintf("Dataset: %s | Before: %d | After: %d | Dropped: %.2f%%",
                  dataset_name, total_rows_before, total_rows_after, dropped_pct))
  return(cleaned_data)
}

# --- Load raw data (same preprocessing as 01_model.R) ---
message("\n--- Loading and preprocessing data ---")

RRST_trial_data <- readr::read_csv(file.path(BASE_DIR, "data", "cleaned", "RRST.csv"),
                                    show_col_types = FALSE) %>%
  filter(ConfResp >= 0 & ConfResp <= 100) %>%
  mutate(
    Conf = ConfResp / 100,
    Stimulus = Stim,
    Drug = factor(drugs,
                  levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    Accuracy = factor(Resp, levels = c(1, 0), labels = c("Correct", "Incorrect"))
  )
RRST_trial_data <- clean_and_report(RRST_trial_data, "RRST_trial_data")

HRD_trl_data <- read.csv(file.path(BASE_DIR, "data", "cleaned", "HRD.csv")) %>%
  mutate(
    Conf = Confidence / 100,
    drugs = as.factor(drugs),
    Drug = factor(drugs,
                  levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    ResponseCorrect = factor(ResponseCorrect,
                             levels = c("True", "False"),
                             labels = c("Correct", "Incorrect")),
    BPM_scaled = scale(listenBPM)
  ) %>%
  filter(subject != "sub_4049")
HRD_trl_data <- clean_and_report(HRD_trl_data, "HRD_trl_data")

# --- Load models ---
load_latest_model <- function(pattern, data_dir) {
  model_files <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  if (length(model_files) == 0) {
    stop("No model files found matching pattern: ", pattern, " in ", data_dir)
  }
  details <- file.info(model_files)
  model_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
  message("Loading model from: ", basename(model_file))
  readRDS(model_file)
}

rrst_model <- load_latest_model("model_rrst_confidence_.*\\.rds", DATA_DIR)
hrd_model  <- load_latest_model("model_hrd_confidence_.*\\.rds", DATA_DIR)

# ============================================================
# 1. Confidence Distribution (Boundary Check)
# ============================================================
message("\n--- 1. Confidence Distribution ---")

p_dist_hrd <- ggplot(HRD_trl_data, aes(x = Conf)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  facet_wrap(~ Drug) +
  labs(title = "HRD: Confidence Distribution by Drug", x = "Confidence", y = "Count") +
  theme_classic(base_size = 12)

p_dist_rrst <- ggplot(RRST_trial_data, aes(x = Conf)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  facet_wrap(~ Drug) +
  labs(title = "RRST: Confidence Distribution by Drug", x = "Confidence", y = "Count") +
  theme_classic(base_size = 12)

p_dist_combined <- p_dist_hrd / p_dist_rrst +
  plot_annotation(title = "Confidence Distributions (Boundary Check for Ordered Beta)")

dist_filename <- sprintf("diag_confidence_distribution_%s.png", TIMESTAMP)
ggsave(file.path(FIGS_DIR, dist_filename), p_dist_combined,
       width = 10, height = 8, dpi = 300, bg = "white")
message("Saved distribution plot to: ", dist_filename)

# Report boundary proportions
for (task_name in c("HRD", "RRST")) {
  d <- if (task_name == "HRD") HRD_trl_data else RRST_trial_data
  n <- nrow(d)
  n_zero <- sum(d$Conf == 0)
  n_one <- sum(d$Conf == 1)
  message(sprintf("%s: N=%d, Conf==0: %d (%.1f%%), Conf==1: %d (%.1f%%)",
                  task_name, n, n_zero, 100 * n_zero / n, n_one, 100 * n_one / n))
}

# ============================================================
# 2. DHARMa Residual Diagnostics
# ============================================================
message("\n--- 2. DHARMa Residual Diagnostics ---")

run_dharma <- function(model, model_name) {
  tryCatch({
    sim_res <- DHARMa::simulateResiduals(model, n = 500)
    outfile <- file.path(FIGS_DIR, sprintf("diag_dharma_%s_%s.png", tolower(model_name), TIMESTAMP))
    png(outfile, width = 10, height = 5, units = "in", res = 300)
    plot(sim_res, title = paste(model_name, "Confidence Model - DHARMa Residuals"))
    dev.off()
    message("Saved DHARMa plot for ", model_name, " to: ", basename(outfile))
  }, error = function(e) {
    message("DHARMa failed for ", model_name, ": ", e$message)
    message("Trying with refit = FALSE...")
    tryCatch({
      sim_res <- DHARMa::simulateResiduals(model, n = 500, refit = FALSE)
      outfile <- file.path(FIGS_DIR, sprintf("diag_dharma_%s_%s.png", tolower(model_name), TIMESTAMP))
      png(outfile, width = 10, height = 5, units = "in", res = 300)
      plot(sim_res, title = paste(model_name, "Confidence Model - DHARMa Residuals (no refit)"))
      dev.off()
      message("Saved DHARMa plot (no refit) for ", model_name, " to: ", basename(outfile))
    }, error = function(e2) {
      message("DHARMa also failed with refit=FALSE for ", model_name, ": ", e2$message)
    })
  })
}

run_dharma(rrst_model, "RRST")
run_dharma(hrd_model, "HRD")

# ============================================================
# 3. VIF (Collinearity) Check
# ============================================================
message("\n--- 3. VIF Check ---")

message("RRST VIF:")
vif_rrst <- performance::check_collinearity(rrst_model)
print(vif_rrst)

message("\nHRD VIF:")
vif_hrd <- performance::check_collinearity(hrd_model)
print(vif_hrd)

# ============================================================
# 4. Random Effects Inspection
# ============================================================
message("\n--- 4. Random Effects ---")

plot_ranefs <- function(model, model_name) {
  re <- ranef(model)$cond$subject
  re$subject <- rownames(re)

  re_long <- re %>%
    tidyr::pivot_longer(cols = -subject, names_to = "term", values_to = "estimate")

  # Flag outliers (> 2 SD from mean within each term)
  re_long <- re_long %>%
    group_by(term) %>%
    mutate(
      mean_est = mean(estimate, na.rm = TRUE),
      sd_est = sd(estimate, na.rm = TRUE),
      is_outlier = abs(estimate - mean_est) > 2 * sd_est
    ) %>%
    ungroup()

  outliers <- re_long %>% filter(is_outlier)
  if (nrow(outliers) > 0) {
    message(model_name, " - Potential random effect outliers (>2 SD):")
    print(outliers %>% select(subject, term, estimate))
  } else {
    message(model_name, " - No random effect outliers detected.")
  }

  p <- ggplot(re_long, aes(x = estimate, y = reorder(subject, estimate))) +
    geom_point(aes(color = is_outlier), size = 1.5) +
    scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "red")) +
    facet_wrap(~ term, scales = "free_x") +
    labs(title = paste(model_name, "- Subject Random Effects"),
         x = "Random Effect Estimate", y = "Subject") +
    theme_classic(base_size = 10) +
    theme(axis.text.y = element_text(size = 5),
          legend.position = "none")

  return(p)
}

p_ranef_rrst <- plot_ranefs(rrst_model, "RRST")
p_ranef_hrd <- plot_ranefs(hrd_model, "HRD")

p_ranef_combined <- p_ranef_rrst / p_ranef_hrd

ranef_filename <- sprintf("diag_subject_ranefs_%s.png", TIMESTAMP)
ggsave(file.path(FIGS_DIR, ranef_filename), p_ranef_combined,
       width = 12, height = 14, dpi = 300, bg = "white")
message("Saved random effects plot to: ", ranef_filename)

message("\nDiagnostics complete.")
