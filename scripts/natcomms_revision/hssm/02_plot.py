#!/usr/bin/env python
"""
Generate plots from saved HSSM model results (no re-fitting needed).

Usage:
    python 02_plot.py --task hrd                   # most recent HRD run
    python 02_plot.py --task both                  # most recent of each
    python 02_plot.py --run-id hrd_20260211_173258 # specific run
    python 02_plot.py --list                       # list available runs
"""

import argparse
import re

import arviz as az
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

from hssm_utils import FIGS_DIR, TASK_NAMES, TASK_CODING, load_fitted, list_fitted


# ============================================================
# Helpers
# ============================================================

def _clean_label(var):
    """Shorten bambi C(drug, Treatment('PLACEBO')) names."""
    s = str(var)
    s = s.replace("C(drug, Treatment('PLACEBO'))", "drug")
    return s


def _is_fixed_effect(var):
    """True for population-level parameters (not random effects)."""
    s = str(var)
    return "participant_id" not in s and "_sigma" not in s


# ============================================================
# Posterior plot (with separate BISO / PROP panels)
# ============================================================

def plot_posteriors(traces, task, filepath):
    """Plot posterior distributions with separate panels per drug level.

    The ArviZ traces store drug contrasts as a 3rd dimension:
        v_C(drug, Treatment('PLACEBO')): shape (chains, draws, 2)
        with coordinate ['BISO', 'PROP'] on the drug dim.

    This function iterates over that dimension to create separate
    panels for each drug level.
    """
    posterior = traces.posterior
    task_name = TASK_NAMES[task]
    coding = TASK_CODING[task]

    # Collect fixed-effect variables (skip random effects and sigmas)
    fixed_vars = [v for v in posterior.data_vars if _is_fixed_effect(v)]

    # Build panel specs: (var_name, dim_index, label, color)
    # dim_index = None for scalar params, int for indexed drug levels
    panels = []

    DRUG_COLORS = {
        "BISO": "#377EB8",   # blue
        "PROP": "#E41A1C",   # red
    }
    INTER_COLORS = {
        "BISO": "#7FCDBB",   # light teal
        "PROP": "#FC9272",   # light salmon
    }

    for var in sorted(fixed_vars, key=str):
        da = posterior[var]
        extra_dims = [d for d in da.dims if d not in ("chain", "draw")]

        if not extra_dims:
            # Scalar parameter
            label = _clean_label(var)
            if "Intercept" in str(var):
                color = {"v": "#666666", "a": "#4DAF4A", "t": "#984EA3"
                         }.get(str(var)[0], "#999999")
            elif "alpha" in str(var) or "stim" in str(var):
                color = "#FF7F00"  # orange
            else:
                color = "#999999"
            panels.append((var, None, None, label, color))
        else:
            # Multi-level dimension (e.g. drug[BISO, PROP])
            dim_name = extra_dims[0]
            levels = da.coords[dim_name].values
            is_interaction = ("alpha" in str(var) or "stim" in str(var))
            cmap = INTER_COLORS if is_interaction else DRUG_COLORS

            for i, level in enumerate(levels):
                label = f"{_clean_label(var)}[{level}]"
                color = cmap.get(str(level), "#999999")
                panels.append((var, dim_name, i, label, color))

    n = len(panels)
    if n == 0:
        print("  No fixed-effect variables found.")
        return

    ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for i, (var, dim_name, idx, label, color) in enumerate(panels):
        ax = axes[i]
        da = posterior[var]

        if idx is not None:
            samples = da.isel({dim_name: idx}).values.flatten()
        else:
            samples = da.values.flatten()

        ax.hist(samples, bins=60, alpha=0.7, color=color, density=True)

        mean_val = np.mean(samples)
        hdi = az.hdi(np.array(samples), hdi_prob=0.95)

        if "Intercept" not in str(var):
            ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
            ax.axvspan(hdi[0], hdi[1], alpha=0.12, color=color)
            p_dir = max(np.mean(samples > 0), np.mean(samples < 0))
            subtitle = (f"M={mean_val:.3f}  95% HDI [{hdi[0]:.3f}, {hdi[1]:.3f}]"
                        f"\nP(direction)={p_dir:.2f}")
        else:
            subtitle = f"M={mean_val:.3f}"

        ax.set_title(f"{label}\n{subtitle}", fontsize=9, pad=6)
        ax.tick_params(labelsize=8)

    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"{task_name} DDM Posteriors ({coding})",
                 fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(filepath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {filepath.name}")


# ============================================================
# Trace plot
# ============================================================

def plot_trace(traces, task, filepath):
    """ArviZ trace plot for convergence diagnostics."""
    try:
        # Only plot fixed effects to keep it readable
        fixed_vars = [v for v in traces.posterior.data_vars
                      if _is_fixed_effect(v)]
        az.plot_trace(traces, var_names=fixed_vars)
        plt.savefig(filepath, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close("all")
        print(f"  Saved: {filepath.name}")
    except Exception as e:
        print(f"  Trace plot skipped: {e}")


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Plot HSSM results")
    parser.add_argument("--task", choices=["hrd", "rrst", "both"])
    parser.add_argument("--run-id", help="Specific run ID")
    parser.add_argument("--list", action="store_true",
                        help="List available fitted models")
    args = parser.parse_args()

    if args.list:
        print("\nAvailable fitted models:")
        for t in ["hrd", "rrst"]:
            runs = list_fitted(t)
            if runs:
                print(f"\n  {t.upper()}:")
                for r in runs:
                    print(f"    {r}")
        return

    if not args.task and not args.run_id:
        parser.error("Provide --task or --run-id (or --list)")

    if args.run_id:
        run_ids = [args.run_id]
    else:
        tasks = ["hrd", "rrst"] if args.task == "both" else [args.task]
        run_ids = []
        for t in tasks:
            available = list_fitted(t)
            if available:
                run_ids.append(available[0])  # most recent
            else:
                print(f"  No fitted models for {t.upper()}")

    for rid in run_ids:
        print(f"\n{'=' * 60}")
        print(f"Plotting: {rid}")
        print("=" * 60)

        result = load_fitted(run_id=rid)
        task = result["config"].get("task", rid.split("_")[0])

        if result["traces"] is None:
            print("  No traces found, skipping.")
            continue

        plot_posteriors(
            result["traces"], task,
            FIGS_DIR / f"{rid}_posteriors.png",
        )
        plot_trace(
            result["traces"], task,
            FIGS_DIR / f"{rid}_trace.png",
        )

    print(f"\n====== Done! Plots in {FIGS_DIR} ======")


if __name__ == "__main__":
    main()
