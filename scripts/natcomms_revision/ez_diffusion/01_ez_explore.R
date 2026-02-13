#!/usr/bin/env Rscript
# ============================================================
# EZ-Diffusion Model: Quick Exploration
# ============================================================
# Decomposes RT + accuracy into drift rate (v), boundary
# separation (a), and non-decision time (Ter) for HRD and
# RRST tasks across drug conditions (PLACEBO, BISO, PROP).
#
# Uses hausekeep::fit_ezddm (Hause Lin) for parameter estimation.
#
# Reference: Wagenmakers, van der Maas, & Grasman (2007).
#   An EZ-diffusion model for response time and accuracy.
#   Psychonomic Bulletin & Review, 14(1), 3-22.
#
# Limitation: EZ assumes no starting-point bias (z = a/2).
#   If participants have a response bias (e.g., tendency to
#   say "More" in HRD), this is absorbed into drift rate.
# ============================================================

if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org")
pacman::p_load(dplyr, readr, tidyr, ggplot2, here, lmerTest, emmeans, data.table)

# Install hausekeep from GitHub if not available
if (!requireNamespace("hausekeep", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("hauselin/hausekeep")
}
library(hausekeep)

# --- Directory setup (matches project conventions) ---
BASE_DIR    <- here::here()
SCRIPT_DIR  <- file.path(BASE_DIR, "scripts", "natcomms_revision", "ez_diffusion")
OUTPUT_DIR  <- file.path(SCRIPT_DIR, "results")
FIGS_DIR    <- file.path(OUTPUT_DIR, "figs")
DATA_DIR    <- file.path(OUTPUT_DIR, "data")
DOCS_DIR    <- file.path(OUTPUT_DIR, "docs")

dir.create(FIGS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Output directory: ", OUTPUT_DIR)
message("Run Timestamp:    ", TIMESTAMP)


# ============================================================
# 1. HEARTBEAT DISCRIMINATION TASK (HRD)
# ============================================================
message("\n====== HRD: Heartbeat Discrimination Task ======")

hrd_raw <- readr::read_csv(
  file.path(BASE_DIR, "data", "HRD_trial_level_data.csv"),
  show_col_types = FALSE
)

# Filter: interoceptive trials, valid decisions, valid RT
hrd <- hrd_raw %>%
  filter(
    Modality == "Intero",
    Decision %in% c("More", "Less"),
    !is.na(DecisionRT)
  ) %>%
  mutate(
    correct = as.integer(ResponseCorrect),
    drugs   = factor(drugs, levels = c("PLACEBO", "BISO", "PROP"))
  ) %>%
  filter(drugs %in% c("PLACEBO", "BISO", "PROP"))

# Subject alignment with cleaned data
cleaned_hrd_subs <- unique(
  read.csv(file.path(BASE_DIR, "data", "cleaned", "HRD.csv"))$subject
)
hrd <- hrd %>% filter(subject %in% cleaned_hrd_subs)

# Outlier removal (matches hrd_logit pipeline)
hrd <- hrd %>% filter(subject != "sub_4071")

# RT trimming: anticipatory (< 200 ms) and lapse (> 5 s)
n_before <- nrow(hrd)
hrd <- hrd %>% filter(DecisionRT >= 0.2, DecisionRT <= 5)
message(sprintf("  RT trimming: dropped %d / %d trials (%.1f%%)",
                n_before - nrow(hrd), n_before,
                100 * (n_before - nrow(hrd)) / n_before))

message(sprintf("  Final: %d subjects, %d trials",
                length(unique(hrd$subject)), nrow(hrd)))

# Fit EZ-DDM using hausekeep
hrd_ez <- fit_ezddm(
  data      = hrd,
  rts       = "DecisionRT",
  responses = "correct",
  id        = "subject",
  group     = "drugs"
)

message(sprintf("  EZ parameters: %d cells returned", nrow(hrd_ez)))


# ============================================================
# 2. RESPIRATORY RESISTANCE TASK (RRST)
# ============================================================
message("\n====== RRST: Respiratory Resistance Sensation Task ======")

rrst_raw <- readr::read_csv(
  file.path(BASE_DIR, "data", "RRST_trial_level_data.csv"),
  show_col_types = FALSE
)

rrst <- rrst_raw %>%
  filter(!is.na(RT)) %>%
  mutate(
    correct = as.integer(trialSuccess),
    drugs   = factor(drugs, levels = c("PLACEBO", "BISO", "PROP"))
  ) %>%
  filter(drugs %in% c("PLACEBO", "BISO", "PROP"))

# Subject alignment with cleaned data
cleaned_rrst_subs <- unique(
  read.csv(file.path(BASE_DIR, "data", "cleaned", "RRST.csv"))$subject
)
rrst <- rrst %>% filter(subject %in% cleaned_rrst_subs)

# RT trimming
n_before <- nrow(rrst)
rrst <- rrst %>% filter(RT >= 0.2, RT <= 5)
message(sprintf("  RT trimming: dropped %d / %d trials (%.1f%%)",
                n_before - nrow(rrst), n_before,
                100 * (n_before - nrow(rrst)) / n_before))

message(sprintf("  Final: %d subjects, %d trials",
                length(unique(rrst$subject)), nrow(rrst)))

# Fit EZ-DDM using hausekeep
rrst_ez <- fit_ezddm(
  data      = rrst,
  rts       = "RT",
  responses = "correct",
  id        = "subject",
  group     = "drugs"
)

message(sprintf("  EZ parameters: %d cells returned", nrow(rrst_ez)))


# ============================================================
# 3. DESCRIPTIVE STATISTICS
# ============================================================
message("\n====== Descriptive Statistics ======")

describe_ez <- function(ez_data, task_name) {
  ez_data %>%
    as.data.frame() %>%
    group_by(drugs) %>%
    summarise(
      n        = n(),
      v_mean   = mean(v, na.rm = TRUE),   v_sd   = sd(v, na.rm = TRUE),
      a_mean   = mean(a, na.rm = TRUE),   a_sd   = sd(a, na.rm = TRUE),
      Ter_mean = mean(t0_Ter, na.rm = TRUE), Ter_sd = sd(t0_Ter, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(task = task_name, .before = 1)
}

desc_hrd  <- describe_ez(hrd_ez, "HRD")
desc_rrst <- describe_ez(rrst_ez, "RRST")
desc_all  <- bind_rows(desc_hrd, desc_rrst)

message("\nHRD:")
print(as.data.frame(desc_hrd), digits = 3)
message("\nRRST:")
print(as.data.frame(desc_rrst), digits = 3)


# ============================================================
# 4. STATISTICAL TESTS — Drug effects on each DDM parameter
# ============================================================
message("\n====== Statistical Tests (LMM: param ~ drugs + (1|subject)) ======")

run_drug_tests <- function(ez_data, task_name) {
  # fit_ezddm returns: a, v, t0_Ter (plus n, n0, n1, behavioural cols)
  param_map <- c(
    v     = "Drift Rate (v)",
    a     = "Boundary Separation (a)",
    t0_Ter = "Non-Decision Time (Ter)"
  )

  # Convert to data.frame for lmer
  df_base <- as.data.frame(ez_data)

  results <- list()

  for (p in names(param_map)) {
    df <- df_base %>% filter(!is.na(.data[[p]]))

    mod <- lmerTest::lmer(
      as.formula(paste(p, "~ drugs + (1 | subject)")),
      data = df
    )

    anova_res <- anova(mod)
    emm       <- emmeans::emmeans(mod, "drugs")
    pairs_res <- pairs(emm, adjust = "holm")
    eff       <- emmeans::eff_size(emm, sigma = sigma(mod), edf = df.residual(mod))

    results[[p]] <- list(
      param    = p,
      label    = param_map[p],
      model    = mod,
      anova    = anova_res,
      emmeans  = emm,
      pairs    = pairs_res,
      eff_size = eff
    )

    # Print to console
    message(sprintf("\n--- %s: %s ---", task_name, param_map[p]))
    f_row <- anova_res[1, ]
    message(sprintf("  F(%d, %.1f) = %.3f, p = %s",
                    f_row$NumDF, f_row$DenDF, f_row$`F value`,
                    format.pval(f_row$`Pr(>F)`, digits = 3)))
    cat("  Marginal means:\n")
    print(summary(emm))
    cat("  Pairwise (Holm-adjusted):\n")
    print(summary(pairs_res))
    cat("  Effect sizes (Cohen's d):\n")
    print(summary(eff))
  }

  results
}

message("\n------ HRD ------")
hrd_results <- run_drug_tests(hrd_ez, "HRD")

message("\n------ RRST ------")
rrst_results <- run_drug_tests(rrst_ez, "RRST")


# ============================================================
# 5. VISUALIZATIONS
# ============================================================
message("\n====== Generating Plots ======")

# Long-format for plotting
to_long <- function(ez_data, task_name) {
  ez_data %>%
    as.data.frame() %>%
    select(subject, drugs, v, a, t0_Ter) %>%
    rename(Ter = t0_Ter) %>%
    pivot_longer(c(v, a, Ter), names_to = "parameter", values_to = "value") %>%
    mutate(task = task_name)
}

plot_df <- bind_rows(
  to_long(hrd_ez, "HRD"),
  to_long(rrst_ez, "RRST")
) %>%
  filter(!is.na(value)) %>%
  mutate(
    parameter = factor(parameter,
      levels = c("v", "a", "Ter"),
      labels = c("Drift Rate (v)", "Boundary Separation (a)", "Non-Decision Time (Ter)")
    ),
    drugs = factor(drugs, levels = c("PLACEBO", "BISO", "PROP"))
  )

drug_colors <- c(PLACEBO = "#4DAF4A", BISO = "#377EB8", PROP = "#E41A1C")

# --- Plot 1: Violin + box + points ---
p1 <- ggplot(plot_df, aes(x = drugs, y = value, fill = drugs)) +
  geom_violin(alpha = 0.35, width = 0.7, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.06, alpha = 0.3, size = 0.8) +
  facet_grid(parameter ~ task, scales = "free_y") +
  scale_fill_manual(values = drug_colors) +
  labs(
    title = "EZ-Diffusion Parameters by Drug Condition",
    subtitle = "Wagenmakers et al. (2007) decomposition of RT + accuracy",
    x = NULL, y = "Parameter Value", fill = "Drug"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(FIGS_DIR, sprintf("ez_violin_%s.png", TIMESTAMP)),
       p1, width = 8, height = 9, dpi = 300, bg = "white")

# --- Plot 2: Within-subject spaghetti (paired lines + group mean +/- SE) ---
p2 <- ggplot(plot_df, aes(x = drugs, y = value)) +
  geom_line(aes(group = subject), alpha = 0.12, color = "grey50") +
  geom_point(aes(color = drugs), alpha = 0.35, size = 1) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               linewidth = 1.2, color = "black") +
  stat_summary(aes(group = 1), fun = mean, geom = "point",
               size = 3, color = "black") +
  stat_summary(aes(group = 1), fun.data = mean_se, geom = "errorbar",
               width = 0.12, linewidth = 0.8, color = "black") +
  facet_grid(parameter ~ task, scales = "free_y") +
  scale_color_manual(values = drug_colors) +
  labs(
    title = "Within-Subject EZ-Diffusion Parameter Changes",
    subtitle = "Grey lines = individual subjects; black = group mean +/- SE",
    x = NULL, y = "Parameter Value", color = "Drug"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(FIGS_DIR, sprintf("ez_spaghetti_%s.png", TIMESTAMP)),
       p2, width = 8, height = 9, dpi = 300, bg = "white")

message("  Plots saved to: ", FIGS_DIR)


# ============================================================
# 6. SAVE DATA
# ============================================================
readr::write_csv(as.data.frame(hrd_ez),
  file.path(DATA_DIR, sprintf("hrd_ez_parameters_%s.csv", TIMESTAMP)))
readr::write_csv(as.data.frame(rrst_ez),
  file.path(DATA_DIR, sprintf("rrst_ez_parameters_%s.csv", TIMESTAMP)))
readr::write_csv(desc_all,
  file.path(DATA_DIR, sprintf("ez_descriptives_%s.csv", TIMESTAMP)))


# ============================================================
# 7. FULL TEXT REPORT
# ============================================================
report_file <- file.path(DOCS_DIR, sprintf("ez_diffusion_report_%s.txt", TIMESTAMP))
sink(report_file)

cat("================================================================\n")
cat("EZ-DIFFUSION MODEL — EXPLORATORY ANALYSIS\n")
cat("Run: ", TIMESTAMP, "\n")
cat("================================================================\n")
cat("\nMethod: EZ-diffusion (Wagenmakers et al., 2007) via hausekeep::fit_ezddm\n")
cat("Scaling: s = 0.1\n")
cat("RT trimming: 0.2 s < RT < 5 s\n")
cat("Statistics: LMM — param ~ drugs + (1 | subject), Holm-adjusted pairwise\n")
cat("\nNote: EZ assumes symmetric starting point (no response bias).\n")
cat("      Parameters are computed per subject x drug cell.\n")

for (task in list(
  list(name = "HRD (Heartbeat Discrimination)", data = hrd, ez = hrd_ez, res = hrd_results),
  list(name = "RRST (Respiratory Resistance)",  data = rrst, ez = rrst_ez, res = rrst_results)
)) {
  cat("\n\n================================================================\n")
  cat(task$name, "\n")
  cat(sprintf("  Subjects: %d  |  Trials: %d\n",
              length(unique(task$data$subject)), nrow(task$data)))
  cat("================================================================\n")

  for (r in task$res) {
    cat(sprintf("\n--- %s ---\n", r$label))
    cat("\nEstimated marginal means:\n")
    print(summary(r$emmeans))
    cat("\nType III ANOVA:\n")
    print(r$anova)
    cat("\nPairwise contrasts (Holm):\n")
    print(summary(r$pairs))
    cat("\nEffect sizes (Cohen's d):\n")
    print(summary(r$eff_size))
  }
}

sink()
message("\n  Report saved to: ", report_file)

message("\n====== Done! ======")
message("All outputs in: ", OUTPUT_DIR)
