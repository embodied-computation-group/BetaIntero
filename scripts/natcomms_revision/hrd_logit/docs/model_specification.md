# Model Specification: Interoceptive Choice Under Beta-Blockade

Reference script: `01_model.R`

---

## 1. Research Question

Does pharmacological manipulation of cardiac signalling (beta-adrenergic blockade) alter interoceptive accuracy, as measured by participants' ability to judge whether an auditory tone rate is faster or slower than their own heartbeat?

The primary test is the **Drug x Listen BPM interaction**: whether the slope relating perceived heart rate (listenBPM) to choice differs across drug conditions.

---

## 2. Data Source

**Input file:** `data/HRD_trial_level_data.csv`

Each row is a single trial from a heartbeat discrimination task. Key columns:

| Column | Description |
|---|---|
| `Modality` | Task modality; only `"Intero"` trials are retained |
| `Decision` | Participant response: `"More"` or `"Less"` |
| `listenBPM` | Heart rate (BPM) the participant listened to during the trial |
| `responseBPM` | Tone rate (BPM) presented for comparison |
| `subject` | Participant ID (e.g., `sub_4071`) |
| `drugs` | Drug condition: `PLACEBO`, `BISO` (bisoprolol), or `PROP` (propranolol) |
| `nTrials` | Trial number within the session |
| `visit` | Visit number (1, 2, or 3) |

---

## 3. Data Preparation

### 3.1 Filtering

1. **Modality filter:** Only interoceptive (`Modality == "Intero"`) trials are retained.
2. **Valid responses:** Trials are kept only if `Decision` is `"More"` or `"Less"`. Trials with missing `listenBPM` or `responseBPM` are dropped.
3. **Drug conditions:** Only the three conditions of interest (`PLACEBO`, `BISO`, `PROP`) are retained.

### 3.2 Outcome Coding

The binary outcome variable `Response_binary` is coded as:
- 1 = `"More"` (participant judged the tones as faster than their heartbeat)
- 0 = `"Less"`

### 3.3 Factor Coding

- **`drugs`**: Dummy-coded factor with `PLACEBO` as the reference level. The two estimated contrasts are BISO vs PLACEBO and PROP vs PLACEBO.
- **`visit`**: Converted to an unordered factor (dummy-coded). Visit 1 is the reference level.

### 3.4 Outlier Removal

Subject `sub_4071` is excluded prior to model fitting. Diagnostics (see `03_diagnostics.R`) identified this participant as having an implausible individual listenBPM slope (> 70 log-odds units), indicating near-perfect separation or data quality issues.

---

## 4. Centering

Continuous BPM predictors are **mean-centered within subject**: each participant's grand mean (pooled across all drug sessions) is subtracted from the trial-level value. No division by SD is performed, so coefficients remain on the original BPM scale.

| Variable | Raw source | Centering | Interpretation |
|---|---|---|---|
| `listenBPM_centered` | `listenBPM` | subject mean subtracted | Trial-level deviation from that participant's overall mean listen BPM |
| `responseBPM_centered` | `responseBPM` | subject mean subtracted | Trial-level deviation from that participant's overall mean response BPM |
| `nTrials` | `nTrials` | none (raw) | Raw trial number within the session |

No between-subject component or separate between-subject regressor is computed.

---

## 5. Model Specification

### 5.1 Estimator

Generalised linear mixed model fitted via maximum likelihood using `glmmTMB`.

### 5.2 Family and Link

Binomial family with logit link. Coefficients are on the log-odds scale.

### 5.3 Fixed Effects

```
Response_binary ~ drugs * listenBPM_centered
                + responseBPM_centered
                + nTrials
                + visit
```

This expands to the following terms:

| Term | Interpretation |
|---|---|
| `(Intercept)` | Log-odds of responding "More" at the mean of all continuous predictors, for PLACEBO at Visit 1 |
| `drugsBISO` | Shift in intercept (log-odds) for bisoprolol vs placebo |
| `drugsPROP` | Shift in intercept (log-odds) for propranolol vs placebo |
| `listenBPM_centered` | Slope of listen BPM on log-odds of "More" (under placebo) |
| `responseBPM_centered` | Slope of response BPM on log-odds of "More" (averaged across drugs; no interaction) |
| `nTrials` | Effect of trial number (learning/fatigue control) |
| `visit` (2 vs 1) | Session-level shift at Visit 2 relative to Visit 1 |
| `visit` (3 vs 1) | Session-level shift at Visit 3 relative to Visit 1 |
| `drugsBISO:listenBPM_centered` | **Key test.** Difference in the listenBPM slope under bisoprolol vs placebo |
| `drugsPROP:listenBPM_centered` | **Key test.** Difference in the listenBPM slope under propranolol vs placebo |

### 5.4 Random Effects

The target (maximal) random effects structure is:

```
(1 + listenBPM_centered + responseBPM_centered | subject)
```

This estimates:
- A random intercept per subject
- A random slope for `listenBPM_centered` per subject
- A random slope for `responseBPM_centered` per subject
- All pairwise correlations among the three random terms

### 5.5 Convergence Strategy

If the maximal model fails to converge (non-zero convergence code or non-positive-definite Hessian), the script falls back through two simpler structures:

1. **Uncorrelated random slopes** `(1 + listenBPM_centered + responseBPM_centered || subject)` -- drops the random-effect correlation parameters.
2. **Random intercept only** `(1 | subject)` -- drops all random slopes.

Convergence is assessed by checking:
- `fit$convergence == 0` (optimiser reports success)
- `sdr$pdHess == TRUE` (Hessian is positive definite, indicating a proper maximum)

---

## 6. Post-hoc Tests

After model fitting, simple slopes of `listenBPM_centered` are estimated separately for each drug condition using `emmeans::emtrends()`. This answers: "Is the listenBPM slope significantly different from zero within each drug?"

Pairwise comparisons of these slopes (via `pairs()`) directly test whether the listenBPM slope differs between any two drug conditions (e.g., BISO vs PLACEBO, PROP vs PLACEBO, BISO vs PROP).

---

## 7. Outputs

| Output | Location | Content |
|---|---|---|
| Model object | `outputs/data/model_combined_drug_interaction_with_trialN_<timestamp>.rds` | Serialised `glmmTMB` model for downstream scripts |
| Model summary | `outputs/docs/model_summary_combined_drug_interaction_with_trialN_<timestamp>.txt` | Full printed summary including N subjects, N trials, formula, fixed effects, random effects |
| Coefficients CSV | `outputs/data/model_coefficients_combined_with_trialN_<timestamp>.csv` | Fixed effects table (estimate, SE, z, p) plus AIC and timestamp |
| Emmeans slopes | `outputs/data/emmeans_listenBPM_slopes_trialN_<timestamp>.csv` | Simple slope of listenBPM per drug condition |
| Pairwise comparisons | `outputs/data/emmeans_listenBPM_slopes_pairwise_trialN_<timestamp>.csv` | Slope difference tests between drug conditions |

---

## 8. Downstream Scripts

| Script | Depends on | Purpose |
|---|---|---|
| `02_plot.R` | Model `.rds` file | Forest plot of fixed effects (odds ratios); marginal effects of each predictor; Drug x Listen BPM interaction plot |
| `03_diagnostics.R` | Model `.rds` file + raw data | Correlation checks, raw-data effect visualisation, individual subject slope distributions, VIF collinearity check |
| `04_apa_table.R` | Model `.rds` file | APA-formatted fixed effects table exported as DOCX, CSV, and PDF |
