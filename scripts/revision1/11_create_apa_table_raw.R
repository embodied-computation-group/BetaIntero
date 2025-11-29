#!/usr/bin/env Rscript
# Create APA-Style Results Table for BetaIntero Revision 1 (Raw Model)
# Extracts fixed effects from the latest raw model and saves as a Word document and PDF.

if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org")
pacman::p_load(dplyr, readr, glmmTMB, broom.mixed, flextable, officer, here, webshot2, htmltools)

# Setup directories
BASE_DIR <- here::here()
SCRIPT_DIR <- file.path(BASE_DIR, "scripts", "revision1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "outputs")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

# Create tables directory if it doesn't exist
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# Find most recent model file matching the pattern
model_pattern <- "model_raw_interaction_.*\\.rds"
model_files <- list.files(DATA_DIR, pattern = model_pattern, full.names = TRUE)

if (length(model_files) == 0) {
  stop("No raw model files found in ", DATA_DIR)
} else {
  # Sort by modification time to get the most recent
  details <- file.info(model_files)
  model_file <- rownames(details)[order(details$mtime, decreasing = TRUE)[1]]
}

message("Loading model from: ", basename(model_file))
mod <- readRDS(model_file)

# Extract fixed effects using broom.mixed
message("Extracting fixed effects...")
tidy_res <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

# Format for APA Style
results_formatted <- tidy_res %>%
  mutate(
    # Rename terms for readability
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
    # Format numbers
    B = sprintf("%.3f", estimate), # 3 decimals for raw coefficients might be better given small slopes
    SE = sprintf("%.3f", std.error),
    CI_95 = sprintf("[%.3f, %.3f]", conf.low, conf.high),
    z = sprintf("%.2f", statistic),
    # Format p-values
    p = case_when(
      p.value < .001 ~ "< .001",
      p.value < .01 ~ sprintf("%.3f", p.value),
      TRUE ~ sprintf("%.2f", p.value)
    )
  ) %>%
  select(Predictor, B, SE, CI_95, z, p)

# Create Flextable
message("Creating APA table...")
ft <- flextable(results_formatted) %>%
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
  
  # Column widths (adjust as needed)
  width(j = "Predictor", width = 2.5) %>%
  width(j = c("B", "SE", "z", "p"), width = 0.8) %>%
  width(j = "CI_95", width = 1.5) %>%
  
  # Caption
  set_caption(
    caption = as_paragraph(
      as_chunk("Table 1", props = fp_text_lite(bold = TRUE)),
      as_chunk("\n"),
      as_chunk("Fixed Effects for Interoception Choice Model (Raw BPM Analysis)")
    )
  )

# Generate Timestamp for filename
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
doc_filename <- sprintf("apa_table_raw_results_%s.docx", TIMESTAMP)
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
csv_filename <- sprintf("apa_table_raw_results_%s.csv", TIMESTAMP)
csv_file <- file.path(TABLES_DIR, csv_filename)
write.csv(results_formatted, csv_file, row.names = FALSE)
message("Saved CSV table to: ", csv_file)

# Save as PDF
html_filename <- sprintf("apa_table_raw_results_%s.html", TIMESTAMP)
html_file <- file.path(TABLES_DIR, html_filename)

# Generate HTML using htmltools
html_val <- flextable::htmltools_value(ft)
htmltools::save_html(html_val, file = html_file)

pdf_filename <- sprintf("apa_table_raw_results_%s.pdf", TIMESTAMP)
pdf_file <- file.path(TABLES_DIR, pdf_filename)

message("Saving as PDF...")
tryCatch({
  webshot2::webshot(url = html_file, file = pdf_file, vwidth = 800, vheight = 600)
  message("Saved APA table to: ", pdf_file)
  # Optional: Remove the temporary HTML file
  # file.remove(html_file)
}, error = function(e) {
  message("Failed to save PDF. Ensure 'webshot2' is installed and a browser (Chrome/Edge) is available.")
  message("Error details: ", e$message)
})

message("Table generation complete!")
