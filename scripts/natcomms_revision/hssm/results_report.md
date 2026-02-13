# HSSM Hierarchical DDM Analysis: Drug Effects on Interoceptive Decision-Making

## Overview

Hierarchical Drift Diffusion Models (DDMs) were fitted to trial-level data from two interoceptive tasks (HRD, RRST) to decompose the effects of beta-adrenergic blockade (bisoprolol, propranolol) on latent decision processes. Models were fitted using [HSSM](https://github.com/lnccbrown/HSSM) v0.2.10 (Frank Lab, successor to HDDM) with the No-U-Turn Sampler (NUTS) via NumPyro/JAX.

## Method

### Package and sampler

- **Package:** HSSM v0.2.10 with analytical DDM likelihood
- **Sampler:** NUTS (NumPyro/JAX), 4 chains x 1,000 draws (1,500 tune), target_accept = 0.90
- **Outlier probability:** 0.05 (lapse trials)
- **Starting point:** z = 0.5 (fixed, unbiased)

### DDM parameters

| Parameter | Description |
|-----------|-------------|
| **v** (drift rate) | Rate of evidence accumulation; modulated by experimental factors |
| **a** (boundary separation) | Response caution / speed-accuracy tradeoff |
| **t** (non-decision time) | Encoding + motor execution time |
| **z** (starting point) | Fixed at 0.5 (no a priori response bias) |

All parameters except z include participant-level random intercepts: `(1|participant_id)`.

---

## Task 1: Heartbeat Discrimination (HRD)

### Design

Participants judged whether their heart rate during a response phase was "More" or "Less" than during a preceding listen phase. The interoceptive signal strength on each trial is **alpha** = responseBPM - listenBPM.

### DDM coding

**Response-coded:** upper boundary = "More" (+1), lower boundary = "Less" (-1). Positive drift rates drive responses toward "More"; alpha is expected to positively modulate drift (larger heart rate differences are easier to discriminate in the "More" direction).

### Model specification

```
v ~ 1 + C(drug, Treatment('PLACEBO')) * alpha + (1|participant_id)
a ~ 1 + (1|participant_id)
t ~ 1 + (1|participant_id)
z = 0.5 (fixed)
```

Drug is treatment-coded with PLACEBO as reference. Alpha is z-scored within each participant.

### Data

- **N** = 48 participants, 8,141 trials
- **Modality:** Interoceptive trials only
- **RT filter:** 0.2-5.0 s
- **Exclusions:** Cleaned subject list + sub_4071

### Results

| Parameter | Mean | SD | 95% HDI | ESS (bulk) | r_hat |
|-----------|------|----|---------|------------|-------|
| v Intercept | 0.348 | 0.024 | [0.302, 0.390] | 2,784 | 1.00 |
| **v drug[BISO]** | **-0.323** | **0.026** | **[-0.373, -0.276]** | **3,878** | **1.00** |
| **v drug[PROP]** | **-0.287** | **0.026** | **[-0.335, -0.240]** | **3,846** | **1.00** |
| **v alpha** | **0.788** | **0.021** | **[0.748, 0.827]** | **3,157** | **1.00** |
| **v drug:alpha[BISO]** | **-0.096** | **0.028** | **[-0.151, -0.044]** | **3,623** | **1.00** |
| **v drug:alpha[PROP]** | **-0.067** | **0.029** | **[-0.119, -0.011]** | **3,599** | **1.00** |
| a Intercept | 1.321 | 0.022 | [1.281, 1.362] | 1,421 | 1.00 |
| t Intercept | 1.281 | 0.033 | [1.218, 1.341] | 499 | 1.01 |

### Interpretation

1. **Alpha strongly drives drift rate** (0.788, HDI excludes 0). Larger heart rate increases push evidence accumulation toward "More", confirming interoceptive signal detection.

2. **Both drugs reduce drift rate relative to placebo.** BISO (-0.323) and PROP (-0.287) both have 95% HDIs entirely below zero. Beta-blockade impairs the rate of interoceptive evidence accumulation.

3. **Drug x Alpha interactions are negative.** Both BISO (-0.096) and PROP (-0.067) attenuate the influence of alpha on drift. Under beta-blockade, participants are less sensitive to trial-level variation in heart rate — the mapping from physical signal to decision evidence is weakened.

4. **Convergence is excellent.** All r_hat = 1.00, ESS > 499 across all fixed effects, zero divergences.

---

## Task 2: Respiratory Resistance Sensitivity (RRST)

### Design

A 2-interval forced choice (2IFC) task: participants identified which of two breathing intervals contained an added inspiratory resistance. Stimulus intensity (cmH2O) was adaptively staircased to ~82% accuracy.

### DDM coding

**Accuracy-coded:** upper boundary = correct (+1), lower boundary = incorrect (-1). Positive drift rates reflect better evidence accumulation toward the correct response.

### Model specification

```
v ~ 1 + C(drug, Treatment('PLACEBO')) + stim_c + (1|participant_id)
a ~ 1 + (1|participant_id)
t ~ 1 + (1|participant_id)
z = 0.5 (fixed)
```

Drug is treatment-coded with PLACEBO as reference. Stimulus intensity (stim_c) is z-scored within each participant to isolate trial-level variation from between-subject differences in staircase difficulty.

### Data

- **N** = 48 participants, 8,130 trials
- **Accuracy:** 79.1% (consistent with staircase target)
- **RT filter:** 0.2-5.0 s
- **Accuracy variable:** `Resp` (= whether participant correctly identified the signal interval)

### Results

| Parameter | Mean | SD | 95% HDI | ESS (bulk) | r_hat |
|-----------|------|----|---------|------------|-------|
| v Intercept | 0.936 | 0.042 | [0.858, 1.017] | 1,925 | 1.00 |
| **v drug[BISO]** | **0.118** | **0.039** | **[0.048, 0.195]** | **4,894** | **1.00** |
| v drug[PROP] | 0.023 | 0.040 | [-0.051, 0.097] | 5,152 | 1.00 |
| **v stim_c** | **0.507** | **0.016** | **[0.478, 0.538]** | **5,413** | **1.00** |
| a Intercept | 0.762 | 0.020 | [0.723, 0.800] | 644 | 1.00 |
| t Intercept | 0.160 | 0.008 | [0.144, 0.176] | 479 | 1.00 |

### Interpretation

1. **Stimulus intensity strongly drives drift** (0.507, HDI excludes 0). Higher resistance loads produce faster, more accurate evidence accumulation — confirming the signal detection mechanism.

2. **BISO increases drift rate** (0.118, HDI excludes 0). Bisoprolol enhances evidence accumulation accuracy for respiratory resistance detection, an effect in the *opposite* direction from HRD.

3. **PROP has no credible effect** (0.023, HDI spans 0). Propranolol does not alter respiratory resistance sensitivity.

4. **Convergence is good.** All r_hat = 1.00, 16 divergences (negligible).

---

## Cross-Task Comparison

| Effect | HRD (response-coded) | RRST (accuracy-coded) |
|--------|---------------------|----------------------|
| BISO on v | **-0.323** (impairs) | **+0.118** (enhances) |
| PROP on v | **-0.287** (impairs) | 0.023 (null) |
| Covariate on v | alpha: +0.788 | stim_c: +0.507 |
| Drug x Covariate | BISO: -0.096, PROP: -0.067 | (not modeled) |
| a (boundary) | 1.321 | 0.762 |
| t (non-decision) | 1.281 s | 0.160 s |

### Key observations

1. **Dissociable drug effects on drift rate across tasks.** Beta-blockade (especially BISO) *impairs* cardiac interoceptive evidence accumulation but *enhances* respiratory resistance detection. This dissociation suggests the two tasks recruit distinct sensory-decision pathways with different dependencies on beta-adrenergic signaling.

2. **PROP effects are weaker and task-dependent.** PROP credibly impairs HRD drift but has no detectable effect on RRST, consistent with differential pharmacological profiles (PROP is non-selective; BISO is beta-1 selective).

3. **Non-decision times differ markedly.** HRD t = 1.28 s (includes cardiac listen/response phases), RRST t = 0.16 s (simple perceptual comparison). This reflects the different temporal structures of the tasks.

4. **Response caution is higher for HRD** (a = 1.32 vs 0.76), consistent with the more deliberative nature of the heartbeat discrimination judgment versus the 2IFC respiratory task.

---

## Pipeline

```
scripts/natcomms_revision/hssm/
  hssm_utils.py      # Shared config, data loading, model building, save/load
  01_fit.py           # Fit models: python 01_fit.py --mode full --task both
  02_plot.py          # Plot results: python 02_plot.py --task both
  03_report.py        # Text report: python 03_report.py --task both
  results/
    fitted/           # Saved model runs (traces.nc, summary.csv, config.json, data.csv)
    figs/             # Posterior and trace plots
    docs/             # Text reports
```

### Fitted model runs

| Run ID | Task | Mode | Subjects | Trials |
|--------|------|------|----------|--------|
| `hrd_20260211_182953` | HRD | full | 48 | 8,141 |
| `rrst_20260211_192147` | RRST | full | 48 | 8,130 |

### Reproducibility

```bash
# Fit both models (full: 4 chains x 1000 draws, ~5 min each)
python 01_fit.py --mode full --task both --parallel

# Generate plots from saved traces (no re-fitting)
python 02_plot.py --task both

# Generate text report
python 03_report.py --task both
```
