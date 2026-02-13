# Nature Communications Revision Analyses

Supplementary analyses for the Nature Communications revision, examining the effects of beta-adrenergic blockade (propranolol, bisoprolol vs placebo) on interoceptive processing.

All scripts use `here::here()` for path resolution and save timestamped outputs to a `results/` subdirectory within each pipeline.

---

## Directory Structure

```
natcomms_revision/
├── hrd_logit/           # HRD perception (interoceptive choice)
│   ├── 01_model.R
│   ├── 02_plot.R
│   ├── 03_diagnostics.R
│   ├── 04_apa_table.R
│   ├── docs/
│   └── results/
└── metacognition/       # Metacognitive confidence (RRST + HRD)
    ├── 01_model.R
    ├── 02_plot.R
    ├── 03_diagnostics.R
    ├── 04_apa_table.R
    ├── 05_rrst_interaction.R
    ├── 06_rrst_interaction_diagnostics.R
    ├── 07_rrst_interaction_apa_table.R
    ├── 01b_model_bpm_within_subject.R
    ├── plot_model_vs_data.R
    ├── test_compare.R
    └── results/
```

---

## `hrd_logit/` — HRD Interoceptive Choice Model

**Question:** Does beta-blockade alter interoceptive accuracy on the Heartbeat Discrimination (HRD) task?

**Model:** Binomial GLMM (logit link) fitted with `glmmTMB`.

```
Response_binary ~ drugs * listenBPM_centered + responseBPM_centered + nTrials + visit
                + (1 + listenBPM_centered + responseBPM_centered | subject)
```

- **Data:** `data/HRD_trial_level_data.csv` (interoceptive trials only)
- **Outcome:** Binary choice (More/Less) indicating whether the participant judged tones as faster than their heartbeat
- **Key test:** Drug x listenBPM interaction — whether the psychometric slope differs across drug conditions
- **Centering:** BPM predictors are subject-mean-centered (no scaling)
- **Outlier:** sub_4071 excluded (implausible slope from separation)
- **Post-hoc:** Simple slopes via `emmeans::emtrends()` + pairwise slope comparisons

| Script | Purpose |
|---|---|
| `01_model.R` | Fit model, save coefficients, emmeans |
| `02_plot.R` | Forest plot, marginal effects, Drug x BPM interaction |
| `03_diagnostics.R` | Correlation checks, subject slopes, VIF |
| `04_apa_table.R` | APA-format fixed effects table (DOCX, CSV, PDF) |

See [`docs/model_specification.md`](hrd_logit/docs/model_specification.md) for full model documentation.

---

## `metacognition/` — Metacognitive Confidence Models

**Question:** Does beta-blockade affect metacognitive confidence on the RRST (respiratory) and HRD (cardiac) tasks?

**Model:** Ordered beta regression (`family = ordbeta()`) fitted with `glmmTMB`, appropriate for confidence ratings bounded at [0, 1].

### Core pipeline (scripts 01–04)

Two models are fit in `01_model.R`:

**RRST:**
```
Conf ~ Stimulus_c + Drug * Accuracy + (1 + Stimulus_c + Accuracy | subject)
```
- Stimulus_c = subject-mean-centered stimulus intensity

**HRD:**
```
Conf ~ Drug * ResponseCorrect + BPM_scaled + (1 + ResponseCorrect + BPM_scaled | subject)
```
- BPM_scaled = within-subject z-scored listenBPM
- sub_4049 excluded

| Script | Purpose |
|---|---|
| `01_model.R` | Fit RRST + HRD confidence models, emmeans, coefficients |
| `02_plot.R` | Predicted confidence plots + pairwise contrast forest plots |
| `03_diagnostics.R` | DHARMa residuals, VIF, boundary checks, random effect outliers |
| `04_apa_table.R` | APA-format tables for both models + HRD post-hoc (DOCX, CSV, PDF) |

### RRST interaction analysis (scripts 05–07)

Extends the RRST model with a three-way Stimulus x Drug x Accuracy interaction to test whether stimulus intensity effects on confidence vary by drug and accuracy.

```
Conf ~ Stimulus_c * Drug * Accuracy + (1 + Stimulus_c + Accuracy | subject)
```

| Script | Purpose |
|---|---|
| `05_rrst_interaction.R` | Fit three-way interaction model, emmeans, emtrends |
| `06_rrst_interaction_diagnostics.R` | DHARMa, VIF, random effects for interaction model |
| `07_rrst_interaction_apa_table.R` | APA tables for interaction fixed effects + slope contrasts |

### Auxiliary scripts

| Script | Purpose |
|---|---|
| `01b_model_bpm_within_subject.R` | Variant of 01 with within-subject z-scored BPM (now matches main 01) |
| `plot_model_vs_data.R` | Overlay model predictions on raw RRST data |
| `test_compare.R` | Validation: compare old (Analysis.Rmd) vs new pipeline output |

---

## Running Order

Each pipeline is independent. Within each, run scripts numerically:

```
# HRD perception
Rscript scripts/natcomms_revision/hrd_logit/01_model.R
Rscript scripts/natcomms_revision/hrd_logit/02_plot.R
Rscript scripts/natcomms_revision/hrd_logit/03_diagnostics.R
Rscript scripts/natcomms_revision/hrd_logit/04_apa_table.R

# Metacognitive confidence
Rscript scripts/natcomms_revision/metacognition/01_model.R
Rscript scripts/natcomms_revision/metacognition/02_plot.R
Rscript scripts/natcomms_revision/metacognition/03_diagnostics.R
Rscript scripts/natcomms_revision/metacognition/04_apa_table.R

# RRST interaction (optional, depends on 01_model.R data)
Rscript scripts/natcomms_revision/metacognition/05_rrst_interaction.R
Rscript scripts/natcomms_revision/metacognition/06_rrst_interaction_diagnostics.R
Rscript scripts/natcomms_revision/metacognition/07_rrst_interaction_apa_table.R
```

Scripts 02–04 depend on the model `.rds` files saved by 01. Each script automatically loads the most recently saved model by timestamp.

## Output Structure

Both pipelines write to `results/` with the same layout:

```
results/
├── data/     # Model .rds objects, coefficient CSVs, emmeans CSVs
├── docs/     # Model summary text files
├── figs/     # PNG plots
└── tables/   # APA tables (DOCX, CSV, HTML, PDF)
```

All output filenames include a `YYYYMMDD_HHMMSS` timestamp. The `results/` directories are gitignored.
