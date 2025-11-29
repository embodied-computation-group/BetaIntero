#!/usr/bin/env Rscript
# Create APA-Style Results Table for Final Raw Model (no 3-way interaction)

if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org")
pacman::p_load(dplyr, broom.mixed, flextable, officer, here, webshot2, htmltools)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "revision1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# Load final model (no 3-way interaction) - pick newest saved model
model_files <- list.files(DATA_DIR, pattern = "^model_raw_interaction_.*\\.rds$", full.names = TRUE)
if (length(model_files) == 0) stop("No model files found in data directory: ", DATA_DIR)
model_file <- model_files[order(file.info(model_files)$mtime, decreasing = TRUE)][1]
message("Loading model from: ", model_file)
mod <- readRDS(model_file)

# Prefer reading the pre-saved coefficients CSV (created by model script) if available
coef_files <- list.files(file.path(OUTPUT_DIR, "docs"), pattern = "^model_coefficients_.*\\.csv$", full.names = TRUE)
if (length(coef_files) > 0) {
  coef_file <- coef_files[order(file.info(coef_files)$mtime, decreasing = TRUE)][1]
  message("Loading coefficients from: ", coef_file)
  coef_tab <- readr::read_csv(coef_file, show_col_types = FALSE)
  if (!"term" %in% names(coef_tab)) {
    coef_tab <- coef_tab %>% tibble::rownames_to_column(var = "term")
  }
  tidy_res <- coef_tab %>%
    dplyr::rename(estimate = Estimate, std.error = `Std. Error`, statistic = `z value`, p.value = `Pr(>|z|)`) %>%
    dplyr::mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error
    )
} else {
  # Fallback: Extract fixed effects from summary (robust for glmmTMB objects)
  message("No pre-saved coefficients CSV found; extracting from model object...")
  mod_sum <- summary(mod)
  if (!is.list(mod_sum) || is.null(mod_sum$coefficients)) stop("Unexpected model summary structure; cannot extract coefficients")
  coef_tab <- as.data.frame(mod_sum$coefficients$cond) %>% tibble::rownames_to_column(var = "term")
  coef_tab <- coef_tab %>%
    dplyr::rename(estimate = Estimate, std.error = `Std. Error`, statistic = `z value`, p.value = `Pr(>|z|)`) %>%
    dplyr::mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error
    )
  tidy_res <- coef_tab
}

# Format for APA Style
results_formatted <- tidy_res %>%
  mutate(
    Predictor = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "drugsBISO" ~ "Drug: Bisoprolol (vs Placebo)",
      term == "drugsPROP" ~ "Drug: Propranolol (vs Placebo)",
      term == "listenBPM_centered" ~ "Listen BPM (Centered)",
      term == "responseBPM_centered" ~ "Response BPM (Centered)",
      term == "nTrials" ~ "Trial Number",
      term == "visitvisit_0002" ~ "Visit 2 (vs Visit 1)",
      term == "visitvisit_0003" ~ "Visit 3 (vs Visit 1)",
      term == "drugsBISO:listenBPM_centered" ~ "Bisoprolol × Listen BPM",
      term == "drugsPROP:listenBPM_centered" ~ "Propranolol × Listen BPM",
      TRUE ~ term
    ),
    B = sprintf("%.3f", estimate),
    SE = sprintf("%.3f", std.error),
    CI_95 = sprintf("[%.3f, %.3f]", conf.low, conf.high),
    z = sprintf("%.2f", statistic),
    p = case_when(
      p.value < .001 ~ "< .001",
      p.value < .01 ~ sprintf("%.3f", p.value),
      TRUE ~ sprintf("%.2f", p.value)
    )
  ) %>%
  select(Predictor, B, SE, CI_95, z, p)

# Create Flextable
ft <- flextable(results_formatted) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 11, part = "all") %>%
  bold(part = "header") %>%
  align(j = c("B", "SE", "CI_95", "z", "p"), align = "center", part = "all") %>%
  align(j = "Predictor", align = "left", part = "all") %>%
  hline_top(border = fp_border(width = 2), part = "all") %>%
  hline_bottom(border = fp_border(width = 2), part = "body") %>%
  hline(i = 1, border = fp_border(width = 1), part = "header") %>%
  width(j = "Predictor", width = 2.5) %>%
  width(j = c("B", "SE", "z", "p"), width = 0.8) %>%
  width(j = "CI_95", width = 1.5) %>%
  set_caption(
    caption = as_paragraph(
      as_chunk("Table 1", props = fp_text_lite(bold = TRUE)),
      as_chunk("\n"),
      as_chunk("Fixed Effects for Interoception Choice Model (Final Raw BPM Model)")
    )
  )

# Generate Timestamp for filename
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
doc_filename <- sprintf("apa_table_final_raw_results_%s.docx", TIMESTAMP)
output_file <- file.path(TABLES_DIR, doc_filename)

# Save to Word Document
doc <- read_docx() %>%
  body_add_par("") %>%
  body_add_flextable(ft) %>%
  body_add_par("") %>%
  body_add_par("Note. B = unstandardized coefficient (log-odds scale); SE = standard error; CI_95 = 95% confidence interval. Predictors are unstandardized (Raw BPM). Reference level for Drug is Placebo. Model: Response ~ drugs * listenBPM + responseBPM + nTrials + visit + (RE).", 
               style = "Normal")

print(doc, target = output_file)
message("Saved APA table to: ", output_file)

# Also save as CSV
csv_filename <- sprintf("apa_table_final_raw_results_%s.csv", TIMESTAMP)
csv_file <- file.path(TABLES_DIR, csv_filename)
write.csv(results_formatted, csv_file, row.names = FALSE)
message("Saved CSV table to: ", csv_file)

# Save as PDF
html_filename <- sprintf("apa_table_final_raw_results_%s.html", TIMESTAMP)
html_file <- file.path(TABLES_DIR, html_filename)
html_val <- flextable::htmltools_value(ft)
htmltools::save_html(html_val, file = html_file)
pdf_filename <- sprintf("apa_table_final_raw_results_%s.pdf", TIMESTAMP)
pdf_file <- file.path(TABLES_DIR, pdf_filename)
tryCatch({
  webshot2::webshot(url = html_file, file = pdf_file, vwidth = 800, vheight = 600)
  message("Saved APA table to: ", pdf_file)
}, error = function(e) {
  message("Failed to save PDF. Ensure 'webshot2' is installed and a browser (Chrome/Edge) is available.")
  message("Error details: ", e$message)
})

message("Table generation complete!")
