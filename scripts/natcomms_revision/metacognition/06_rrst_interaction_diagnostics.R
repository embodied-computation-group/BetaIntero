#!/usr/bin/env Rscript
# Diagnostics for RRST Stimulus * Drug * Accuracy interaction model
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

# --- Load interaction model (most recent) ---
model_files <- list.files(DATA_DIR,
                          pattern = "model_rrst_confidence_interaction_.*\\.rds",
                          full.names = TRUE)
if (length(model_files) == 0) stop("No interaction model files found in ", DATA_DIR)
details <- file.info(model_files)
model_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
message("Loading model: ", basename(model_file))
mod <- readRDS(model_file)

# --- Recreate data (same preprocessing as 05_rrst_interaction.R) ---
RRST_trial_data <- readr::read_csv(
  file.path(BASE_DIR, "data", "cleaned", "RRST.csv"),
  show_col_types = FALSE
) %>%
  filter(ConfResp >= 0 & ConfResp <= 100) %>%
  mutate(
    Conf = ConfResp / 100,
    Stimulus = Stim,
    Drug = factor(drugs,
                  levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    Accuracy = factor(Resp, levels = c(1, 0), labels = c("Correct", "Incorrect"))
  ) %>%
  drop_na() %>%
  group_by(subject) %>%
  mutate(Stimulus_c = Stimulus - mean(Stimulus, na.rm = TRUE)) %>%
  ungroup()

message("Data: ", nrow(RRST_trial_data), " trials, ",
        length(unique(RRST_trial_data$subject)), " subjects")

# ============================================================
# 1. Confidence Distribution (Boundary Check)
# ============================================================
message("\n--- 1. Confidence Distribution ---")

n <- nrow(RRST_trial_data)
n_zero <- sum(RRST_trial_data$Conf == 0)
n_one <- sum(RRST_trial_data$Conf == 1)
message(sprintf("RRST: N=%d, Conf==0: %d (%.1f%%), Conf==1: %d (%.1f%%)",
                n, n_zero, 100 * n_zero / n, n_one, 100 * n_one / n))

p_dist <- ggplot(RRST_trial_data, aes(x = Conf)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  facet_grid(Accuracy ~ Drug) +
  labs(title = "RRST Interaction Model: Confidence Distribution",
       subtitle = "By Drug and Accuracy (boundary check for ordered beta)",
       x = "Confidence", y = "Count") +
  theme_classic(base_size = 12)

ggsave(file.path(FIGS_DIR,
                 sprintf("diag_interaction_conf_dist_%s.png", TIMESTAMP)),
       p_dist, width = 10, height = 6, dpi = 300, bg = "white")
message("Saved distribution plot.")

# ============================================================
# 2. DHARMa Residual Diagnostics
# ============================================================
message("\n--- 2. DHARMa Residual Diagnostics ---")

tryCatch({
  sim_res <- DHARMa::simulateResiduals(mod, n = 500)
  outfile <- file.path(FIGS_DIR,
                       sprintf("diag_interaction_dharma_%s.png", TIMESTAMP))
  png(outfile, width = 10, height = 5, units = "in", res = 300)
  plot(sim_res, title = "RRST Interaction Model - DHARMa Residuals")
  dev.off()
  message("Saved DHARMa plot to: ", basename(outfile))

  # Additional DHARMa tests
  message("\nDHARMa uniformity test:")
  print(DHARMa::testUniformity(sim_res))

  message("\nDHARMa dispersion test:")
  print(DHARMa::testDispersion(sim_res))

  message("\nDHARMa outlier test:")
  print(DHARMa::testOutliers(sim_res))

}, error = function(e) {
  message("DHARMa failed: ", e$message)
  message("Trying with refit = FALSE...")
  tryCatch({
    sim_res <- DHARMa::simulateResiduals(mod, n = 500, refit = FALSE)
    outfile <- file.path(FIGS_DIR,
                         sprintf("diag_interaction_dharma_%s.png", TIMESTAMP))
    png(outfile, width = 10, height = 5, units = "in", res = 300)
    plot(sim_res, title = "RRST Interaction Model - DHARMa Residuals (no refit)")
    dev.off()
    message("Saved DHARMa plot (no refit) to: ", basename(outfile))
  }, error = function(e2) {
    message("DHARMa also failed with refit=FALSE: ", e2$message)
  })
})

# ============================================================
# 3. VIF (Collinearity) Check
# ============================================================
message("\n--- 3. VIF Check ---")

tryCatch({
  vif_res <- performance::check_collinearity(mod)
  message("VIF results:")
  print(vif_res)
}, error = function(e) {
  message("VIF check failed: ", e$message)
})

# ============================================================
# 4. Random Effects Inspection
# ============================================================
message("\n--- 4. Random Effects ---")

re <- ranef(mod)$cond$subject
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
  message("Potential random effect outliers (>2 SD):")
  print(as.data.frame(outliers %>% select(subject, term, estimate)))
} else {
  message("No random effect outliers detected.")
}

p_ranef <- ggplot(re_long, aes(x = estimate, y = reorder(subject, estimate))) +
  geom_point(aes(color = is_outlier), size = 1.5) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "red")) +
  facet_wrap(~ term, scales = "free_x") +
  labs(title = "RRST Interaction Model - Subject Random Effects",
       x = "Random Effect Estimate", y = "Subject") +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 5),
        legend.position = "none")

ggsave(file.path(FIGS_DIR,
                 sprintf("diag_interaction_ranefs_%s.png", TIMESTAMP)),
       p_ranef, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved random effects plot.")

# ============================================================
# 5. Model Summary Recap
# ============================================================
message("\n--- 5. Model Summary ---")
message("Formula: ", deparse(formula(mod)))
message("AIC: ", AIC(mod))
message("Convergence: ", mod$sdr$pdHess)

message("\nDiagnostics complete.")
