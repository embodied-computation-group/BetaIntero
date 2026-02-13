#!/usr/bin/env Rscript
# APA Table for RRST Stimulus * Drug * Accuracy interaction model
# Fixed effects + emtrends slope contrasts

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, glmmTMB, broom.mixed, emmeans, flextable, officer, here, webshot2)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision", "metacognition")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Table Run Timestamp: ", TIMESTAMP)

# --- Load most recent interaction model ---
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

mod <- load_latest_model("model_rrst_confidence_interaction_.*\\.rds", DATA_DIR)

# --- Recreate data so emmeans/emtrends can access it ---
RE_OUTLIER_SUBJECTS <- c("sub_3004", "sub_3005", "sub_4032", "sub_4054",
                         "sub_4068", "sub_4076", "sub_4079")

RRST_trial_data <- readr::read_csv(
  file.path(BASE_DIR, "data", "cleaned", "RRST.csv"),
  show_col_types = FALSE
) %>%
  dplyr::filter(ConfResp >= 0 & ConfResp <= 100) %>%
  dplyr::mutate(
    Conf = ConfResp / 100,
    Stimulus = Stim,
    Drug = factor(drugs,
                  levels = c("PLACEBO", "PROP", "BISO"),
                  labels = c("placebo", "propranolol", "bisoprolol")),
    Drug = relevel(Drug, ref = "placebo"),
    Accuracy = factor(Resp, levels = c(1, 0), labels = c("Correct", "Incorrect"))
  ) %>%
  tidyr::drop_na() %>%
  dplyr::group_by(subject) %>%
  dplyr::mutate(Stimulus_c = Stimulus - mean(Stimulus, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!subject %in% RE_OUTLIER_SUBJECTS) %>%
  dplyr::filter(Stimulus >= 8)

message("Data: ", nrow(RRST_trial_data), " trials, ",
        length(unique(RRST_trial_data$subject)), " subjects",
        " (Stimulus >= 8)")

# --- Helper: format model to APA table ---
format_apa_table <- function(model, predictor_labels) {
  tidy_res <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

  unmatched <- setdiff(tidy_res$term, names(predictor_labels))
  if (length(unmatched) > 0) {
    message("Warning: unmatched terms (using raw names): ", paste(unmatched, collapse = ", "))
  }

  results_formatted <- tidy_res %>%
    mutate(
      Predictor = ifelse(term %in% names(predictor_labels),
                         predictor_labels[term],
                         term),
      B  = sprintf("%.2f", estimate),
      SE = sprintf("%.2f", std.error),
      CI_95 = sprintf("[%.2f, %.2f]", conf.low, conf.high),
      z  = sprintf("%.2f", statistic),
      p  = case_when(
        p.value < .001 ~ "< .001",
        p.value < .01  ~ sprintf("%.3f", p.value),
        TRUE           ~ sprintf("%.2f", p.value)
      )
    ) %>%
    select(Predictor, B, SE, CI_95, z, p)

  return(results_formatted)
}

# --- Helper: create and save flextable ---
save_apa_flextable <- function(results_df, caption_text, note_text, file_prefix) {
  ft <- flextable(results_df) %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 11, part = "all") %>%
    bold(part = "header") %>%
    align(j = setdiff(names(results_df), "Predictor"), align = "center", part = "all") %>%
    align(j = "Predictor", align = "left", part = "all") %>%
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1), part = "header") %>%
    width(j = "Predictor", width = 3.0) %>%
    set_caption(
      caption = as_paragraph(
        as_chunk(caption_text, props = fp_text_lite(bold = TRUE))
      )
    )

  # Set widths for standard columns if present
  std_cols <- intersect(c("B", "SE", "z", "p"), names(results_df))
  if (length(std_cols) > 0) ft <- width(ft, j = std_cols, width = 0.8)
  if ("CI_95" %in% names(results_df)) ft <- width(ft, j = "CI_95", width = 1.2)

  # Save DOCX
  doc_file <- file.path(TABLES_DIR, sprintf("%s_%s.docx", file_prefix, TIMESTAMP))
  doc <- read_docx() %>%
    body_add_par("") %>%
    body_add_flextable(ft) %>%
    body_add_par("") %>%
    body_add_par(note_text, style = "Normal")
  print(doc, target = doc_file)
  message("Saved DOCX to: ", doc_file)

  # Save CSV
  csv_file <- file.path(TABLES_DIR, sprintf("%s_%s.csv", file_prefix, TIMESTAMP))
  write.csv(results_df, csv_file, row.names = FALSE)
  message("Saved CSV to: ", csv_file)

  # Save PDF via HTML + webshot2
  if (!require("htmltools")) install.packages("htmltools")
  html_file <- file.path(TABLES_DIR, sprintf("%s_%s.html", file_prefix, TIMESTAMP))
  html_val <- flextable::htmltools_value(ft)
  htmltools::save_html(html_val, file = html_file)

  pdf_file <- file.path(TABLES_DIR, sprintf("%s_%s.pdf", file_prefix, TIMESTAMP))
  tryCatch({
    webshot2::webshot(url = html_file, file = pdf_file, vwidth = 800, vheight = 600)
    message("Saved PDF to: ", pdf_file)
  }, error = function(e) {
    message("Failed to save PDF: ", e$message)
  })
}

# ============================================================
# Table 1: Fixed Effects
# ============================================================
message("\n--- Fixed Effects Table ---")

interaction_labels <- c(
  "(Intercept)"                                    = "Intercept",
  "Stimulus_c"                                     = "Stimulus Intensity (centered)",
  "Drugpropranolol"                                = "Drug: Propranolol",
  "Drugbisoprolol"                                 = "Drug: Bisoprolol",
  "AccuracyIncorrect"                              = "Accuracy: Incorrect",
  "Stimulus_c:Drugpropranolol"                     = "Stimulus \u00d7 Propranolol",
  "Stimulus_c:Drugbisoprolol"                      = "Stimulus \u00d7 Bisoprolol",
  "Stimulus_c:AccuracyIncorrect"                   = "Stimulus \u00d7 Incorrect",
  "Drugpropranolol:AccuracyIncorrect"              = "Propranolol \u00d7 Incorrect",
  "Drugbisoprolol:AccuracyIncorrect"               = "Bisoprolol \u00d7 Incorrect",
  "Stimulus_c:Drugpropranolol:AccuracyIncorrect"   = "Stimulus \u00d7 Propranolol \u00d7 Incorrect",
  "Stimulus_c:Drugbisoprolol:AccuracyIncorrect"    = "Stimulus \u00d7 Bisoprolol \u00d7 Incorrect"
)

fixef_table <- format_apa_table(mod, interaction_labels)

save_apa_flextable(
  fixef_table,
  caption_text = "Fixed Effects for RRST Confidence: Stimulus \u00d7 Drug \u00d7 Accuracy Interaction (Ordered Beta Regression)",
  note_text = paste0(
    "Note. B = unstandardized coefficient (logit link); SE = standard error; ",
    "CI_95 = 95% Wald confidence interval. Model: Conf ~ Stimulus_c * Drug * Accuracy + ",
    "(1 + Stimulus_c * Accuracy || subject). Family: ordered beta regression. ",
    "Reference levels: Drug = Placebo, Accuracy = Correct. ",
    "Stimulus_c is subject-mean-centered stimulus intensity. ",
    "7 RE-outlier subjects excluded (sub_3004, sub_3005, sub_4032, sub_4054, sub_4068, sub_4076, sub_4079). ",
    "Trials with Stimulus < 8 excluded (sparse low-intensity levels). ",
    "N = 41 subjects, 7309 trials."
  ),
  file_prefix = "apa_table_rrst_interaction"
)

# ============================================================
# Table 2: Emtrends (Stimulus slopes by Drug x Accuracy)
# ============================================================
message("\n--- Emtrends Table ---")

stim_trends <- emtrends(mod, ~ Drug | Accuracy, var = "Stimulus_c",
                        data = RRST_trial_data)
trend_summary <- as.data.frame(summary(stim_trends))

trend_formatted <- trend_summary %>%
  mutate(
    Predictor = paste0(Drug, " (", Accuracy, ")"),
    Slope = sprintf("%.3f", Stimulus_c.trend),
    SE = sprintf("%.3f", SE),
    CI_95 = sprintf("[%.3f, %.3f]", asymp.LCL, asymp.UCL)
  ) %>%
  select(Predictor, Slope, SE, CI_95)

save_apa_flextable(
  trend_formatted,
  caption_text = "Stimulus-Confidence Slopes (Emtrends) by Drug and Accuracy",
  note_text = "Note. Slopes represent the linear effect of subject-centered stimulus intensity on confidence (logit scale) at each Drug \u00d7 Accuracy combination. 95% asymptotic Wald confidence intervals.",
  file_prefix = "apa_table_rrst_interaction_emtrends"
)

# ============================================================
# Table 3: Pairwise contrasts of slopes
# ============================================================
message("\n--- Pairwise Slope Contrasts Table ---")

stim_trend_pairs <- pairs(stim_trends)
pairs_summary <- as.data.frame(summary(stim_trend_pairs))

pairs_formatted <- pairs_summary %>%
  mutate(
    Contrast = as.character(contrast),
    Level = as.character(Accuracy),
    Estimate = sprintf("%.3f", estimate),
    SE = sprintf("%.3f", SE),
    z = sprintf("%.2f", z.ratio),
    p = case_when(
      p.value < .001 ~ "< .001",
      p.value < .01  ~ sprintf("%.3f", p.value),
      TRUE           ~ sprintf("%.2f", p.value)
    )
  ) %>%
  select(Contrast, Level, Estimate, SE, z, p)

# Use a simpler save for this table (no "Predictor" column)
ft_pairs <- flextable(pairs_formatted) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 11, part = "all") %>%
  bold(part = "header") %>%
  align(j = c("Estimate", "SE", "z", "p"), align = "center", part = "all") %>%
  align(j = c("Contrast", "Level"), align = "left", part = "all") %>%
  hline_top(border = fp_border(width = 2), part = "all") %>%
  hline_bottom(border = fp_border(width = 2), part = "body") %>%
  hline(i = 1, border = fp_border(width = 1), part = "header") %>%
  width(j = "Contrast", width = 2.2) %>%
  width(j = "Level", width = 1.0) %>%
  width(j = c("Estimate", "SE", "z", "p"), width = 0.8) %>%
  set_caption(
    caption = as_paragraph(
      as_chunk("Pairwise Contrasts of Stimulus-Confidence Slopes by Drug (Tukey-Adjusted)",
               props = fp_text_lite(bold = TRUE))
    )
  )

pairs_doc_file <- file.path(TABLES_DIR, sprintf("apa_table_rrst_interaction_slope_contrasts_%s.docx", TIMESTAMP))
doc_p <- read_docx() %>%
  body_add_par("") %>%
  body_add_flextable(ft_pairs) %>%
  body_add_par("") %>%
  body_add_par("Note. Pairwise comparisons of stimulus-confidence slopes between drug conditions within each accuracy level. P-values adjusted using the Tukey method.",
               style = "Normal")
print(doc_p, target = pairs_doc_file)
message("Saved slope contrasts DOCX to: ", pairs_doc_file)

pairs_csv_file <- file.path(TABLES_DIR, sprintf("apa_table_rrst_interaction_slope_contrasts_%s.csv", TIMESTAMP))
write.csv(pairs_formatted, pairs_csv_file, row.names = FALSE)
message("Saved slope contrasts CSV to: ", pairs_csv_file)

message("\nTable generation complete!")
