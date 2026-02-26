"""
Shared utilities for the HSSM DDM analysis pipeline.

Config, data loading, model building, fitting, and persistence.
"""

import json
import warnings
from datetime import datetime
from pathlib import Path

import cloudpickle

import arviz as az
import hssm
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================
# Paths
# ============================================================
SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parents[2]  # c:\Users\Micah\Desktop\BetaIntero
RESULTS_DIR = SCRIPT_DIR / "results"
FITTED_DIR = RESULTS_DIR / "fitted"
FIGS_DIR = RESULTS_DIR / "figs"
DOCS_DIR = RESULTS_DIR / "docs"

for _d in [FITTED_DIR, FIGS_DIR, DOCS_DIR]:
    _d.mkdir(parents=True, exist_ok=True)

# ============================================================
# Sampler configuration
# ============================================================
MODE_SETTINGS = {
    "debug": dict(sampler="nuts_numpyro", chains=1, draws=100, tune=100,
                  cores=1, target_accept=0.8),
    "quick": dict(sampler="nuts_numpyro", chains=4, draws=500, tune=500,
                  cores=4, target_accept=0.85),
    "full":  dict(sampler="nuts_numpyro", chains=4, draws=1000, tune=1500,
                  cores=4, target_accept=0.9),
}
P_OUTLIER = 0.05

# ============================================================
# Data loading
# ============================================================

def load_hrd_data():
    """Load and prepare HRD data (response-coded DDM).

    Returns DataFrame with columns:
        rt, response, participant_id, drug, alpha
    """
    print("\n====== Loading HRD Data (response-coded) ======")

    hrd_raw = pd.read_csv(BASE_DIR / "data" / "HRD_trial_level_data.csv")
    hrd = hrd_raw[
        (hrd_raw["Modality"] == "Intero")
        & (hrd_raw["Decision"].isin(["More", "Less"]))
        & (hrd_raw["DecisionRT"].notna())
        & (hrd_raw["drugs"].isin(["PLACEBO", "BISO", "PROP"]))
    ].copy()

    cleaned_subs = pd.read_csv(
        BASE_DIR / "data" / "cleaned" / "HRD.csv"
    )["subject"].unique()
    hrd = hrd[hrd["subject"].isin(cleaned_subs)]
    hrd = hrd[hrd["subject"] != "sub_4071"]

    hrd = hrd[(hrd["DecisionRT"] >= 0.2) & (hrd["DecisionRT"] <= 5.0)]

    hrd["rt"] = hrd["DecisionRT"].astype(float)
    hrd["response"] = np.where(hrd["Decision"] == "More", 1.0, -1.0)
    hrd["participant_id"] = hrd["subject"]
    hrd["drug"] = pd.Categorical(
        hrd["drugs"], categories=["PLACEBO", "BISO", "PROP"]
    )
    hrd["alpha_raw"] = hrd["Alpha"].astype(float)
    hrd["alpha"] = hrd.groupby("subject")["alpha_raw"].transform(
        lambda x: (x - x.mean()) / x.std()
    )

    data = hrd[
        ["rt", "response", "participant_id", "drug", "alpha"]
    ].dropna().reset_index(drop=True)

    print(f"  Subjects: {data['participant_id'].nunique()}")
    print(f"  Trials:   {len(data)}")
    print(f"  P(More):  {(data['response'] == 1).mean():.3f}")
    print(f"  Alpha range: [{data['alpha'].min():.2f}, "
          f"{data['alpha'].max():.2f}]")
    return data


def load_rrst_data():
    """Load and prepare RRST data (accuracy-coded DDM).

    Returns DataFrame with columns:
        rt, response, participant_id, drug, stim_c
    """
    print("\n====== Loading RRST Data (accuracy-coded) ======")

    rrst_raw = pd.read_csv(BASE_DIR / "data" / "RRST_trial_level_data.csv")
    rrst = rrst_raw[
        (rrst_raw["RT"].notna())
        & (rrst_raw["drugs"].isin(["PLACEBO", "BISO", "PROP"]))
    ].copy()

    cleaned_subs = pd.read_csv(
        BASE_DIR / "data" / "cleaned" / "RRST.csv"
    )["subject"].unique()
    rrst = rrst[rrst["subject"].isin(cleaned_subs)]

    rrst = rrst[(rrst["RT"] >= 0.2) & (rrst["RT"] <= 5.0)]

    rrst["rt"] = rrst["RT"].astype(float)
    rrst["response"] = np.where(rrst["Resp"] == 1, 1.0, -1.0)
    rrst["participant_id"] = rrst["subject"]
    rrst["drug"] = pd.Categorical(
        rrst["drugs"], categories=["PLACEBO", "BISO", "PROP"]
    )
    rrst["stim_c"] = rrst.groupby("subject")["Stim"].transform(
        lambda x: (x - x.mean()) / x.std()
    )

    data = rrst[
        ["rt", "response", "participant_id", "drug", "stim_c"]
    ].dropna().reset_index(drop=True)

    print(f"  Subjects: {data['participant_id'].nunique()}")
    print(f"  Trials:   {len(data)}")
    print(f"  Accuracy: {(data['response'] == 1).mean():.3f}")
    return data


TASK_LOADERS = {"hrd": load_hrd_data, "rrst": load_rrst_data}

# ============================================================
# Model building
# ============================================================
TASK_FORMULAS = {
    "hrd": {
        "v": "v ~ 1 + C(drug, Treatment('PLACEBO')) * alpha + (1|participant_id)",
        "a": "a ~ 1 + (1|participant_id)",
        "t": "t ~ 1 + (1|participant_id)",
        "z": "0.5 (fixed)",
    },
    "rrst": {
        "v": "v ~ 1 + C(drug, Treatment('PLACEBO')) + stim_c + (1|participant_id)",
        "a": "a ~ 1 + (1|participant_id)",
        "t": "t ~ 1 + (1|participant_id)",
        "z": "0.5 (fixed)",
    },
}
TASK_CODING = {"hrd": "response-coded", "rrst": "accuracy-coded"}
TASK_NAMES = {"hrd": "HRD", "rrst": "RRST"}
TASK_CATEGORICALS = {
    "hrd": {"drug": ["PLACEBO", "BISO", "PROP"]},
    "rrst": {"drug": ["PLACEBO", "BISO", "PROP"]},
}


def build_model(task, data):
    """Build HSSM model for the given task.

    Args:
        task: 'hrd' or 'rrst'
        data: prepared DataFrame from load_*_data()

    Returns:
        hssm.HSSM model (not yet fitted)
    """
    formulas = TASK_FORMULAS[task]

    include = [
        {"name": "v", "formula": formulas["v"], "link": "identity"},
        {"name": "a", "formula": formulas["a"]},
        {"name": "t", "formula": formulas["t"]},
    ]

    model = hssm.HSSM(
        data=data, model="ddm", p_outlier=P_OUTLIER,
        include=include, z=0.5,
    )
    print("Model initialized successfully.")
    return model


# ============================================================
# Fitting
# ============================================================

def fit_model(model, cfg):
    """Fit HSSM model with numpyro NUTS.

    Handles bambi/xarray version incompatibility in log_likelihood
    post-processing (the traces are still valid).

    Args:
        model: HSSM model object
        cfg: dict from MODE_SETTINGS

    Returns:
        ArviZ InferenceData (traces)
    """
    print(f"  Running {cfg['sampler']} ({cfg['chains']} chains x "
          f"{cfg['draws']} draws, tune={cfg['tune']})...")
    try:
        return model.sample(
            sampler=cfg["sampler"], cores=cfg["cores"],
            chains=cfg["chains"], draws=cfg["draws"],
            tune=cfg["tune"], target_accept=cfg["target_accept"],
        )
    except (ValueError, Exception) as e:
        if "different number of dimensions" in str(e):
            print("  NOTE: log-likelihood post-processing failed "
                  "(bambi/xarray compat). Traces are still valid.")
            return getattr(model, "_inference_obj", None)
        raise


# ============================================================
# Persistence: save / load fitted results
# ============================================================

def save_fitted(model, data, task, mode, cfg, timestamp=None):
    """Save fitted model results to disk.

    Creates:  results/fitted/{task}_{timestamp}/
        traces.nc   - ArviZ InferenceData (all posterior samples)
        summary.csv - Parameter summary (mean, sd, HDI, ESS, r_hat)
        config.json - Sampler config + model formulas
        data.csv    - Prepared data used for fitting
        model.pkl   - Cloudpickled HSSM model (for posterior predictive plotting)

    Returns:
        run_id (str): e.g. 'hrd_20260211_180000'
    """
    ts = timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    run_id = f"{task}_{ts}"
    outdir = FITTED_DIR / run_id
    outdir.mkdir(parents=True, exist_ok=True)

    # Traces
    traces = model.traces if hasattr(model, "traces") and model.traces is not None else None
    if traces is not None:
        az.to_netcdf(traces, str(outdir / "traces.nc"))

    # Summary
    try:
        summary = model.summary()
        summary.to_csv(outdir / "summary.csv")
    except Exception as e:
        print(f"  Warning: Could not generate summary: {e}")

    # Config
    config = {
        "task": task,
        "mode": mode,
        "timestamp": ts,
        "run_id": run_id,
        "sampler": cfg,
        "formulas": TASK_FORMULAS[task],
        "coding": TASK_CODING[task],
        "n_subjects": int(data["participant_id"].nunique()),
        "n_trials": len(data),
    }
    with open(outdir / "config.json", "w") as f:
        json.dump(config, f, indent=2)

    # Data
    data.to_csv(outdir / "data.csv", index=False)

    # Model object (for posterior predictive plotting)
    _inf = getattr(model, "_inference_obj", None)
    try:
        model._inference_obj = None  # strip traces (already saved as .nc)
        with open(outdir / "model.pkl", "wb") as f:
            cloudpickle.dump(model, f)
        print(f"  Model object saved.")
    except Exception as e:
        print(f"  Warning: Could not save model object: {e}")
    finally:
        model._inference_obj = _inf  # always restore

    print(f"\n  Saved to: {outdir}")
    return run_id


def load_fitted(task=None, run_id=None):
    """Load fitted model results from disk.

    Args:
        task: load most recent run for this task (e.g. 'hrd')
        run_id: load a specific run (e.g. 'hrd_20260211_180000')

    Returns:
        dict with keys: traces, summary, config, data, model, run_id, outdir
    """
    if run_id:
        outdir = FITTED_DIR / run_id
    elif task:
        dirs = sorted(FITTED_DIR.glob(f"{task}_*"),
                       key=lambda p: p.name, reverse=True)
        if not dirs:
            raise FileNotFoundError(f"No fitted models for task '{task}' "
                                    f"in {FITTED_DIR}")
        outdir = dirs[0]
        run_id = outdir.name
    else:
        raise ValueError("Provide either task or run_id")

    if not outdir.exists():
        raise FileNotFoundError(f"Not found: {outdir}")

    print(f"  Loading: {outdir.name}")

    result = {"run_id": run_id, "outdir": outdir}

    # Traces
    nc = outdir / "traces.nc"
    result["traces"] = az.from_netcdf(str(nc)) if nc.exists() else None

    # Summary
    csv = outdir / "summary.csv"
    result["summary"] = pd.read_csv(csv, index_col=0) if csv.exists() else None

    # Config
    cfg = outdir / "config.json"
    if cfg.exists():
        with open(cfg) as f:
            result["config"] = json.load(f)
    else:
        result["config"] = {}

    # Data
    dcsv = outdir / "data.csv"
    result["data"] = pd.read_csv(dcsv) if dcsv.exists() else None

    # Re-encode categoricals lost in CSV round-trip
    task_name = result["config"].get("task", run_id.split("_")[0])
    if result["data"] is not None and task_name in TASK_CATEGORICALS:
        for col, cats in TASK_CATEGORICALS[task_name].items():
            if col in result["data"].columns:
                result["data"][col] = pd.Categorical(
                    result["data"][col], categories=cats
                )

    # Model object (for posterior predictive / model cartoon plotting)
    result["model"] = None
    pkl = outdir / "model.pkl"
    if pkl.exists():
        try:
            with open(pkl, "rb") as f:
                result["model"] = cloudpickle.load(f)
            result["model"].restore_traces(result["traces"])
            print(f"  Model object loaded from model.pkl")
        except Exception as e:
            print(f"  Warning: Could not load model.pkl: {e}")
            result["model"] = None

    # Fallback: rebuild model from data + config
    if result["model"] is None and result["data"] is not None and result["traces"] is not None:
        try:
            print(f"  Rebuilding model from data + config...")
            result["model"] = build_model(task_name, result["data"])
            result["model"].restore_traces(result["traces"])
            print(f"  Model rebuilt and traces restored.")
        except Exception as e:
            print(f"  Warning: Could not rebuild model: {e}")
            result["model"] = None

    return result


def list_fitted(task=None):
    """List available fitted model runs.

    Returns list of run_id strings, most recent first.
    """
    pattern = f"{task}_*" if task else "*"
    dirs = sorted(FITTED_DIR.glob(pattern),
                   key=lambda p: p.name, reverse=True)
    return [d.name for d in dirs if d.is_dir()]
