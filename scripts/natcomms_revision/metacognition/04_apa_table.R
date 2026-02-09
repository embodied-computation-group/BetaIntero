#!/usr/bin/env Rscript
# Create APA-Style Results Tables for Metacognitive Confidence Models (RRST + HRD)
# Extracts fixed effects and saves as Word document, CSV, and PDF.

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, glmmTMB, broom.mixed, flextable, officer, here, webshot2)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "natcomms_revision", "metacognition")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

# Create tables directory if it doesn't exist
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# Generate Timestamp
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
message("Table Run Timestamp: ", TIMESTAMP)

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

# --- Helper: format model to APA table ---
format_apa_table <- function(model, predictor_labels) {
  tidy_res <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

  # Log any terms not in the label mapping
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
    # Font settings
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 11, part = "all") %>%
    # Header formatting
    bold(part = "header") %>%
    # Alignment
    align(j = c("B", "SE", "CI_95", "z", "p"), align = "center", part = "all") %>%
    align(j = "Predictor", align = "left", part = "all") %>%
    # Borders (APA style)
    hline_top(border = fp_border(width = 2), part = "all") %>%
    hline_bottom(border = fp_border(width = 2), part = "body") %>%
    hline(i = 1, border = fp_border(width = 1), part = "header") %>%
    # Column widths
    width(j = "Predictor", width = 2.5) %>%
    width(j = c("B", "SE", "z", "p"), width = 0.8) %>%
    width(j = "CI_95", width = 1.2) %>%
    # Caption
    set_caption(
      caption = as_paragraph(
        as_chunk(caption_text, props = fp_text_lite(bold = TRUE))
      )
    )

  # Save to Word Document
  doc_file <- file.path(TABLES_DIR, sprintf("%s_%s.docx", file_prefix, TIMESTAMP))
  doc <- read_docx() %>%
    body_add_par("") %>%
    body_add_flextable(ft) %>%
    body_add_par("") %>%
    body_add_par(note_text, style = "Normal")
  print(doc, target = doc_file)
  message("Saved DOCX to: ", doc_file)

  # Save as CSV
  csv_file <- file.path(TABLES_DIR, sprintf("%s_%s.csv", file_prefix, TIMESTAMP))
  write.csv(results_df, csv_file, row.names = FALSE)
  message("Saved CSV to: ", csv_file)

  # Save as PDF via HTML + webshot2
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
# RRST Table
# ============================================================
message("\n--- RRST APA Table ---")

rrst_labels <- c(
  "(Intercept)"                          = "Intercept",
  "Stimulus"                             = "Stimulus Intensity",
  "Drugpropranolol"                      = "Drug: Propranolol (vs Placebo)",
  "Drugbisoprolol"                       = "Drug: Bisoprolol (vs Placebo)",
  "AccuracyIncorrect"                    = "Accuracy: Incorrect (vs Correct)",
  "Drugpropranolol:AccuracyIncorrect"    = "Propranolol \u00d7 Incorrect",
  "Drugbisoprolol:AccuracyIncorrect"     = "Bisoprolol \u00d7 Incorrect"
)

rrst_table <- format_apa_table(rrst_model, rrst_labels)

save_apa_flextable(
  rrst_table,
  caption_text = "Fixed Effects for RRST Metacognitive Confidence Model (Ordered Beta Regression)",
  note_text = "Note. B = unstandardized coefficient (mu link, log-odds scale); SE = standard error; CI_95 = 95% confidence interval. Model family: ordered beta regression (ordbeta). Reference level for Drug is Placebo; reference level for Accuracy is Correct.",
  file_prefix = "apa_table_rrst_confidence"
)

# ============================================================
# HRD Table
# ============================================================
message("\n--- HRD APA Table ---")

hrd_labels <- c(
  "(Intercept)"                                    = "Intercept",
  "Drugpropranolol"                                = "Drug: Propranolol (vs Placebo)",
  "Drugbisoprolol"                                 = "Drug: Bisoprolol (vs Placebo)",
  "ResponseCorrectIncorrect"                       = "Accuracy: Incorrect (vs Correct)",
  "BPM_scaled"                                     = "Listen BPM (Scaled)",
  "Drugpropranolol:ResponseCorrectIncorrect"       = "Propranolol \u00d7 Incorrect",
  "Drugbisoprolol:ResponseCorrectIncorrect"        = "Bisoprolol \u00d7 Incorrect"
)

hrd_table <- format_apa_table(hrd_model, hrd_labels)

save_apa_flextable(
  hrd_table,
  caption_text = "Fixed Effects for HRD Metacognitive Confidence Model (Ordered Beta Regression)",
  note_text = "Note. B = unstandardized coefficient (mu link, log-odds scale); SE = standard error; CI_95 = 95% confidence interval. Model family: ordered beta regression (ordbeta). Reference level for Drug is Placebo; reference level for ResponseCorrect is Correct. Subject sub_4049 excluded.",
  file_prefix = "apa_table_hrd_confidence"
)

# ============================================================
# HRD Post-hoc Table (Pairwise Comparisons)
# ============================================================
message("\n--- HRD Post-hoc Table ---")

# Load pairwise results from 01_model.R output
pairs_pattern <- "emmeans_hrd_pairwise_by_accuracy_.*\\.csv"
pairs_files <- list.files(DATA_DIR, pattern = pairs_pattern, full.names = TRUE)

if (length(pairs_files) > 0) {
  details <- file.info(pairs_files)
  pairs_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
  message("Loading pairwise results from: ", basename(pairs_file))
  pairs_df <- readr::read_csv(pairs_file, show_col_types = FALSE)

  # Format for APA
  pairs_formatted <- pairs_df %>%
    mutate(
      Contrast = contrast,
      Level = ResponseCorrect,
      Estimate = sprintf("%.2f", estimate),
      SE = sprintf("%.2f", SE),
      z = sprintf("%.2f", z.ratio),
      p = case_when(
        p.value < .001 ~ "< .001",
        p.value < .01  ~ sprintf("%.3f", p.value),
        TRUE           ~ sprintf("%.2f", p.value)
      )
    ) %>%
    select(Contrast, Level, Estimate, SE, z, p)

  ft_posthoc <- flextable(pairs_formatted) %>%
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
        as_chunk("Pairwise Comparisons of Drug Effects on HRD Confidence (Holm-Adjusted)",
                 props = fp_text_lite(bold = TRUE))
      )
    )

  # Save posthoc DOCX
  posthoc_doc_file <- file.path(TABLES_DIR, sprintf("apa_table_hrd_posthoc_%s.docx", TIMESTAMP))
  doc_ph <- read_docx() %>%
    body_add_par("") %>%
    body_add_flextable(ft_posthoc) %>%
    body_add_par("") %>%
    body_add_par("Note. Pairwise comparisons of Drug effects within each ResponseCorrect level. P-values adjusted using the Holm method.",
                 style = "Normal")
  print(doc_ph, target = posthoc_doc_file)
  message("Saved posthoc DOCX to: ", posthoc_doc_file)

  # Save posthoc CSV
  posthoc_csv_file <- file.path(TABLES_DIR, sprintf("apa_table_hrd_posthoc_%s.csv", TIMESTAMP))
  write.csv(pairs_formatted, posthoc_csv_file, row.names = FALSE)
  message("Saved posthoc CSV to: ", posthoc_csv_file)

} else {
  message("No pairwise results found. Run 01_model.R first.")
}

message("\nTable generation complete!")
