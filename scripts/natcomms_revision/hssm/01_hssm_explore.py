#!/usr/bin/env python
"""
HSSM Exploratory Analysis - Drug Effects on DDM Parameters
===========================================================
Fits hierarchical drift diffusion models to HRD and RRST data
using HSSM (successor to HDDM, Frank Lab).

HRD:  Response-coded (More=1, Less=-1).
      Drift rate driven by stimulus (Alpha = responseBPM - listenBPM).
      Drug x Alpha interaction tests whether beta-blockers modulate
      the mapping from interoceptive signal to evidence accumulation.

RRST: Accuracy-coded (correct=1, incorrect=-1).
      Drug effect on drift rate tests whether beta-blockers alter
      respiratory evidence quality.

Install (conda recommended):
    conda create -n hssm python=3.11
    conda activate hssm
    conda install -c conda-forge pymc
    pip install hssm

Reference:
    Fengler, Govindarajan, Chen, & Frank (2022).
    https://github.com/lnccbrown/HSSM
"""

import argparse
import subprocess
import sys
import warnings
from datetime import datetime
from pathlib import Path

import arviz as az
import hssm
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")
warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================
# CLI
# ============================================================
parser = argparse.ArgumentParser(description="HSSM DDM analysis")
parser.add_argument(
    "--mode", choices=["debug", "quick", "full"], default="debug",
    help="debug  = numpyro NUTS (1 chain x 100 draws, ~minutes)\n"
         "quick  = numpyro NUTS (4 chains x 500 draws, ~tens of min)\n"
         "full   = numpyro NUTS (4 chains x 1000 draws, ~hours)",
)
parser.add_argument(
    "--task", choices=["hrd", "rrst", "both"], default="both",
    help="Which task to fit (default: both)",
)
parser.add_argument(
    "--parallel", action="store_true",
    help="Fit HRD and RRST in parallel subprocesses (uses ~8 cores)",
)
args = parser.parse_args()

# ============================================================
# Parallel dispatch: spawn two child processes, one per task
# ============================================================
if args.parallel and args.task == "both":
    script = str(Path(__file__).resolve())
    procs = []
    for task in ("hrd", "rrst"):
        cmd = [sys.executable, script,
               "--mode", args.mode, "--task", task]
        print(f"  Spawning: {' '.join(cmd)}")
        procs.append(subprocess.Popen(
            cmd, stdout=sys.stdout, stderr=sys.stderr,
        ))
    codes = [p.wait() for p in procs]
    if any(c != 0 for c in codes):
        print(f"  WARNING: child processes exited with codes {codes}")
    else:
        print("\n  Both tasks completed successfully.")
    sys.exit(max(codes))

# ============================================================
# Configuration (set by --mode)
# ============================================================
MODE_SETTINGS = {
    "debug": dict(sampler="nuts_numpyro", chains=1, draws=100, tune=100,
                  cores=1, target_accept=0.8),
    "quick": dict(sampler="nuts_numpyro", chains=4, draws=500, tune=500,
                  cores=4, target_accept=0.85),
    "full":  dict(sampler="nuts_numpyro", chains=4, draws=1000, tune=1500,
                  cores=4, target_accept=0.9),
}
CFG = MODE_SETTINGS[args.mode]
P_OUTLIER = 0.05

print(f"Mode: {args.mode}  |  Task: {args.task}")
print(f"Sampler: {CFG['sampler']}  |  {CFG['chains']} chains x "
      f"{CFG['draws']} draws (tune={CFG['tune']})")


def fit_model(model):
    """Dispatch sampling based on --mode. Uses numpyro (JAX) for speed."""
    print(f"  Running {CFG['sampler']} ({CFG['chains']} chains x "
          f"{CFG['draws']} draws, tune={CFG['tune']})...")
    try:
        return model.sample(
            sampler=CFG["sampler"],
            cores=CFG["cores"],
            chains=CFG["chains"],
            draws=CFG["draws"],
            tune=CFG["tune"],
            target_accept=CFG["target_accept"],
        )
    except (ValueError, Exception) as e:
        # HSSM auto-computes log_likelihood after sampling, which can
        # fail on certain bambi/xarray version combos. The traces are
        # still stored in model._inference_obj.
        if "different number of dimensions" in str(e):
            print(f"  NOTE: log-likelihood post-processing failed "
                  f"(bambi/xarray compat issue). Traces are still valid.")
            return getattr(model, "_inference_obj", None)
        raise

# ============================================================
# Setup directories
# ============================================================
BASE_DIR = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "results"
FIGS_DIR = OUTPUT_DIR / "figs"
DATA_DIR = OUTPUT_DIR / "data"
DOCS_DIR = OUTPUT_DIR / "docs"

for d in [FIGS_DIR, DATA_DIR, DOCS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Run Timestamp:    {TIMESTAMP}")


# ============================================================
# 1. Load and prepare HRD data (RESPONSE-CODED)
# ============================================================
print("\n====== Loading HRD Data (response-coded) ======")

hrd_raw = pd.read_csv(BASE_DIR / "data" / "HRD_trial_level_data.csv")

hrd = hrd_raw[
    (hrd_raw["Modality"] == "Intero")
    & (hrd_raw["Decision"].isin(["More", "Less"]))
    & (hrd_raw["DecisionRT"].notna())
    & (hrd_raw["drugs"].isin(["PLACEBO", "BISO", "PROP"]))
].copy()

# Subject alignment with cleaned data
cleaned_hrd_subs = pd.read_csv(
    BASE_DIR / "data" / "cleaned" / "HRD.csv"
)["subject"].unique()
hrd = hrd[hrd["subject"].isin(cleaned_hrd_subs)]
hrd = hrd[hrd["subject"] != "sub_4071"]  # known outlier

# RT trimming
hrd = hrd[(hrd["DecisionRT"] >= 0.2) & (hrd["DecisionRT"] <= 5.0)]

# Response coding: More = 1 (upper boundary), Less = -1 (lower)
hrd["rt"] = hrd["DecisionRT"].astype(float)
hrd["response"] = np.where(hrd["Decision"] == "More", 1.0, -1.0)
hrd["participant_id"] = hrd["subject"]

# Drug as categorical with PLACEBO as reference level
hrd["drug"] = pd.Categorical(
    hrd["drugs"], categories=["PLACEBO", "BISO", "PROP"]
)

# Stimulus: Alpha = responseBPM - listenBPM
# Positive Alpha -> tone faster than heartbeat -> evidence for "More"
# Z-score within subject for numerical stability
hrd["alpha_raw"] = hrd["Alpha"].astype(float)
hrd["alpha"] = hrd.groupby("subject")["alpha_raw"].transform(
    lambda x: (x - x.mean()) / x.std()
)

hrd_data = hrd[
    ["rt", "response", "participant_id", "drug", "alpha"]
].dropna().reset_index(drop=True)

print(f"  Subjects: {hrd_data['participant_id'].nunique()}")
print(f"  Trials:   {len(hrd_data)}")
print(f"  P(More):  {(hrd_data['response'] == 1).mean():.3f}")
print(f"  Alpha range: [{hrd_data['alpha'].min():.2f}, "
      f"{hrd_data['alpha'].max():.2f}]")

hrd_data.to_csv(DATA_DIR / f"hrd_prepared_{TIMESTAMP}.csv", index=False)


# ============================================================
# 2. Load and prepare RRST data (ACCURACY-CODED)
# ============================================================
print("\n====== Loading RRST Data (accuracy-coded) ======")

rrst_raw = pd.read_csv(BASE_DIR / "data" / "RRST_trial_level_data.csv")

rrst = rrst_raw[
    (rrst_raw["RT"].notna())
    & (rrst_raw["drugs"].isin(["PLACEBO", "BISO", "PROP"]))
].copy()

# Subject alignment
cleaned_rrst_subs = pd.read_csv(
    BASE_DIR / "data" / "cleaned" / "RRST.csv"
)["subject"].unique()
rrst = rrst[rrst["subject"].isin(cleaned_rrst_subs)]

# RT trimming
rrst = rrst[(rrst["RT"] >= 0.2) & (rrst["RT"] <= 5.0)]

# Accuracy coding: correct = 1, incorrect = -1
rrst["rt"] = rrst["RT"].astype(float)
rrst["response"] = np.where(rrst["trialSuccess"] == 1, 1.0, -1.0)
rrst["participant_id"] = rrst["subject"]
rrst["drug"] = pd.Categorical(
    rrst["drugs"], categories=["PLACEBO", "BISO", "PROP"]
)

# Stimulus intensity (z-scored)
rrst["stim_c"] = (
    (rrst["Stim"] - rrst["Stim"].mean()) / rrst["Stim"].std()
)

rrst_data = rrst[
    ["rt", "response", "participant_id", "drug", "stim_c"]
].dropna().reset_index(drop=True)

print(f"  Subjects: {rrst_data['participant_id'].nunique()}")
print(f"  Trials:   {len(rrst_data)}")
print(f"  Accuracy: {(rrst_data['response'] == 1).mean():.3f}")

rrst_data.to_csv(DATA_DIR / f"rrst_prepared_{TIMESTAMP}.csv", index=False)


# ============================================================
# Helper: save results after fitting a model
# ============================================================
def save_results(model, task_tag, task_name, coding):
    """Save summary, traces, and plots for a fitted model."""
    has_traces = hasattr(model, "traces") and model.traces is not None

    # Summary (works for both Laplace and MCMC)
    try:
        summary = model.summary()
        print(f"\n{task_name} Model Summary:")
        print(summary)
        summary.to_csv(DATA_DIR / f"{task_tag}_summary_{TIMESTAMP}.csv")
    except Exception as e:
        print(f"  Could not generate summary: {e}")
        summary = None

    # Traces + trace plots
    if has_traces:
        az.to_netcdf(
            model.traces,
            str(DATA_DIR / f"{task_tag}_traces_{TIMESTAMP}.nc"),
        )
        try:
            model.plot_trace()
            plt.savefig(
                FIGS_DIR / f"{task_tag}_trace_{TIMESTAMP}.png",
                dpi=150, bbox_inches="tight", facecolor="white",
            )
            plt.close("all")
            print(f"  Trace plot saved.")
        except Exception as e:
            print(f"  Trace plot skipped: {e}")

        # Posterior histograms
        try:
            plot_posteriors(
                model, task_name, coding,
                FIGS_DIR / f"{task_tag}_posteriors_{TIMESTAMP}.png",
            )
            print(f"  Posterior plot saved.")
        except Exception as e:
            print(f"  Posterior plot skipped: {e}")

    return summary


def _clean_varname(var):
    """Shorten bambi C(drug, Treatment('PLACEBO')) names for plot titles."""
    s = str(var)
    # "v_C(drug, Treatment('PLACEBO'))[BISO]" -> "v_drug[BISO]"
    s = s.replace("C(drug, Treatment('PLACEBO'))", "drug")
    # "v_drug:alpha[BISO]" stays as is (readable)
    return s


def plot_posteriors(model, task_name, coding, filepath):
    """Posterior distributions for drug-related and intercept params."""
    posterior = model.traces.posterior
    all_vars = list(posterior.data_vars)

    # ---- Categorise variables ----
    # Main drug effects on v (contain "drug" but NOT "alpha"/"stim")
    v_drug_main = [v for v in all_vars
                   if "drug" in str(v) and "v_" in str(v)
                   and "alpha" not in str(v) and "stim" not in str(v)]
    # Drug x covariate interactions on v
    v_drug_inter = [v for v in all_vars
                    if "drug" in str(v) and "v_" in str(v)
                    and ("alpha" in str(v) or "stim" in str(v))]
    # Main alpha / stim effects on v (no drug)
    v_covariate = [v for v in all_vars
                   if ("alpha" in str(v) or "stim" in str(v))
                   and "v_" in str(v) and "drug" not in str(v)]

    # Ordered panels: intercept, main drug, main covariate, interactions, a, t
    panel_specs = []  # list of (var_name, clean_label, color)

    # v Intercept
    if "v_Intercept" in posterior:
        panel_specs.append(("v_Intercept", "v Intercept", "#666666"))

    # Main drug effects
    for v in sorted(v_drug_main, key=str):
        c = "#377EB8" if "BISO" in str(v) else "#E41A1C"
        panel_specs.append((v, _clean_varname(v), c))

    # Main covariate effects (alpha / stim)
    for v in sorted(v_covariate, key=str):
        panel_specs.append((v, _clean_varname(v), "#FF7F00"))

    # Drug x covariate interactions
    for v in sorted(v_drug_inter, key=str):
        # Light blue for BISO interaction, light red for PROP
        c = "#7FCDBB" if "BISO" in str(v) else "#FC9272"
        panel_specs.append((v, _clean_varname(v), c))

    # a and t intercepts
    if "a_Intercept" in posterior:
        panel_specs.append(("a_Intercept", "a Intercept", "#4DAF4A"))
    if "t_Intercept" in posterior:
        panel_specs.append(("t_Intercept", "t Intercept", "#984EA3"))

    n = len(panel_specs)
    if n == 0:
        print("  No variables found for posterior plot.")
        return

    # Layout: 2 rows if many panels
    ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 3.5 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for i, (var, label, color) in enumerate(panel_specs):
        ax = axes[i]
        s = posterior[var].values.flatten()
        ax.hist(s, bins=60, alpha=0.7, color=color, density=True)

        mean_val = np.mean(s)
        hdi = az.hdi(np.array(s), hdi_prob=0.95)

        # Reference line at zero for effect parameters
        if var not in ("v_Intercept", "a_Intercept", "t_Intercept"):
            ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
            ax.axvspan(hdi[0], hdi[1], alpha=0.12, color=color)
            p_dir = max(np.mean(s > 0), np.mean(s < 0))
            subtitle = (f"M={mean_val:.3f}  95% HDI [{hdi[0]:.3f}, {hdi[1]:.3f}]"
                        f"\nP(direction)={p_dir:.2f}")
        else:
            subtitle = f"M={mean_val:.3f}"

        ax.set_title(f"{label}\n{subtitle}", fontsize=9, pad=6)
        ax.tick_params(labelsize=8)

    # Hide unused axes
    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        f"{task_name} DDM Posteriors ({coding})",
        fontsize=13, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    plt.savefig(filepath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


# ============================================================
# 3. Fit HRD Model
# ============================================================
# Response-coded DDM. Drift rate depends on stimulus (alpha)
# with drug interaction.
#
# v ~ drug * alpha + (1 | participant_id)
#   alpha slope   = perceptual sensitivity
#   drug:alpha    = do drugs modulate this sensitivity?
# a ~ 1 + (1 | participant_id)
# t ~ 1 + (1 | participant_id)
# z = 0.5 (fixed, no response bias)
# ============================================================
hrd_model = None
hrd_summary = None

if args.task in ("hrd", "both"):
    print("\n" + "=" * 60)
    print("Fitting HRD Model")
    print("=" * 60)
    print("  v ~ C(drug, Treatment('PLACEBO')) * alpha + (1|participant_id)")
    print("  a ~ 1 + (1|participant_id)")
    print("  t ~ 1 + (1|participant_id)")
    print("  z = 0.5 (fixed)\n")

    hrd_model = hssm.HSSM(
        data=hrd_data,
        model="ddm",
        p_outlier=P_OUTLIER,
        include=[
            {
                "name": "v",
                "formula": "v ~ 1 + C(drug, Treatment('PLACEBO')) * alpha + (1|participant_id)",
                "link": "identity",
            },
            {
                "name": "a",
                "formula": "a ~ 1 + (1|participant_id)",
            },
            {
                "name": "t",
                "formula": "t ~ 1 + (1|participant_id)",
            },
        ],
        z=0.5,
    )

    fit_model(hrd_model)
    hrd_summary = save_results(
        hrd_model, "hrd", "HRD", "response-coded",
    )


# ============================================================
# 4. Fit RRST Model
# ============================================================
# Accuracy-coded DDM. Positive drift = evidence toward correct.
#
# v ~ drug + stim_c + (1 | participant_id)
# a ~ 1 + (1 | participant_id)
# t ~ 1 + (1 | participant_id)
# z = 0.5 (fixed)
#
# Note: RRST accuracy is near ceiling for many subjects.
# ============================================================
rrst_model = None
rrst_summary = None

if args.task in ("rrst", "both"):
    print("\n" + "=" * 60)
    print("Fitting RRST Model")
    print("=" * 60)
    print("  v ~ C(drug, Treatment('PLACEBO')) + stim_c + (1|participant_id)")
    print("  a ~ 1 + (1|participant_id)")
    print("  t ~ 1 + (1|participant_id)")
    print("  z = 0.5 (fixed)\n")

    rrst_model = hssm.HSSM(
        data=rrst_data,
        model="ddm",
        p_outlier=P_OUTLIER,
        include=[
            {
                "name": "v",
                "formula": "v ~ 1 + C(drug, Treatment('PLACEBO')) + stim_c + (1|participant_id)",
                "link": "identity",
            },
            {
                "name": "a",
                "formula": "a ~ 1 + (1|participant_id)",
            },
            {
                "name": "t",
                "formula": "t ~ 1 + (1|participant_id)",
            },
        ],
        z=0.5,
    )

    fit_model(rrst_model)
    rrst_summary = save_results(
        rrst_model, "rrst", "RRST", "accuracy-coded",
    )


# ============================================================
# 5. Text report
# ============================================================
report_path = DOCS_DIR / f"hssm_report_{TIMESTAMP}.txt"
with open(report_path, "w") as f:
    f.write("=" * 64 + "\n")
    f.write("HSSM - HIERARCHICAL DDM EXPLORATORY ANALYSIS\n")
    f.write(f"Run: {TIMESTAMP}  |  Mode: {args.mode}\n")
    f.write("=" * 64 + "\n\n")

    f.write("Method\n")
    f.write("-" * 40 + "\n")
    f.write("Package: HSSM (Frank Lab, successor to HDDM)\n")
    f.write("Model: DDM with analytical likelihood\n")
    f.write(f"Sampler: {CFG['sampler']}\n")
    f.write(f"Sampling: {CFG['chains']} chains x {CFG['draws']} draws ")
    f.write(f"(tune={CFG['tune']}, "
            f"target_accept={CFG['target_accept']})\n")
    f.write(f"Outlier probability: {P_OUTLIER}\n")
    f.write("Starting point: z = 0.5 (fixed)\n\n")

    f.write("HRD: Response-coded (More=1, Less=-1)\n")
    f.write("  v ~ drug * alpha + (1|participant_id)\n")
    f.write("  alpha = (responseBPM - listenBPM), ")
    f.write("z-scored within subject\n\n")

    f.write("RRST: Accuracy-coded (correct=1, incorrect=-1)\n")
    f.write("  v ~ drug + stim_c + (1|participant_id)\n\n")

    fitted = []
    if hrd_model is not None:
        fitted.append(("HRD (Heartbeat Discrimination)",
                        hrd_model, hrd_summary, hrd_data))
    if rrst_model is not None:
        fitted.append(("RRST (Respiratory Resistance)",
                        rrst_model, rrst_summary, rrst_data))

    for name, model, summary, data in fitted:
        f.write("\n" + "=" * 64 + "\n")
        f.write(f"{name}\n")
        f.write(f"  Subjects: {data['participant_id'].nunique()}")
        f.write(f"  |  Trials: {len(data)}\n")
        f.write("=" * 64 + "\n\n")
        if summary is not None:
            f.write(summary.to_string())
        else:
            f.write("(summary unavailable)\n")
        f.write("\n\n")

print(f"\n  Report saved to: {report_path}")
print(f"\n====== Done! ======")
print(f"All outputs in: {OUTPUT_DIR}")
