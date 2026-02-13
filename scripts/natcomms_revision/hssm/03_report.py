#!/usr/bin/env python
"""
Generate text report from saved HSSM model results.

Usage:
    python 03_report.py --task both                   # most recent
    python 03_report.py --task hrd                    # single task
    python 03_report.py --run-id hrd_20260211_173258  # specific run
"""

import argparse
from datetime import datetime
from pathlib import Path

from hssm_utils import DOCS_DIR, load_fitted, list_fitted


def generate_report(results, output_path):
    """Write text report for one or more fitted model results."""
    with open(output_path, "w") as f:
        f.write("=" * 64 + "\n")
        f.write("HSSM - HIERARCHICAL DDM EXPLORATORY ANALYSIS\n")
        f.write(f"Report: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 64 + "\n\n")

        f.write("Method\n")
        f.write("-" * 40 + "\n")
        f.write("Package: HSSM (Frank Lab, successor to HDDM)\n")
        f.write("Model: DDM with analytical likelihood\n")
        f.write("Starting point: z = 0.5 (fixed)\n\n")

        for r in results:
            cfg = r["config"]
            task = cfg.get("task", "unknown")
            task_label = {"hrd": "HRD (Heartbeat Discrimination)",
                          "rrst": "RRST (Respiratory Resistance)"
                          }.get(task, task.upper())

            f.write("\n" + "=" * 64 + "\n")
            f.write(f"{task_label}\n")
            f.write(f"Run: {r['run_id']}  |  Mode: {cfg.get('mode', '?')}\n")
            f.write("-" * 64 + "\n\n")

            # Sampler
            sc = cfg.get("sampler", {})
            f.write(f"Sampler: {sc.get('sampler', '?')}\n")
            f.write(f"Sampling: {sc.get('chains', '?')} chains x "
                    f"{sc.get('draws', '?')} draws "
                    f"(tune={sc.get('tune', '?')}, "
                    f"target_accept={sc.get('target_accept', '?')})\n")
            f.write(f"Outlier probability: 0.05\n\n")

            # Formulas
            formulas = cfg.get("formulas", {})
            f.write("Model specification:\n")
            for param in ("v", "a", "t", "z"):
                if param in formulas:
                    f.write(f"  {param}: {formulas[param]}\n")
            f.write("\n")

            # Data info
            f.write(f"Subjects: {cfg.get('n_subjects', '?')}")
            f.write(f"  |  Trials: {cfg.get('n_trials', '?')}\n\n")

            # Coding
            f.write(f"Coding: {cfg.get('coding', '?')}\n")
            if task == "hrd":
                f.write("  More = 1 (upper boundary), Less = -1 (lower)\n")
                f.write("  alpha = (responseBPM - listenBPM), "
                        "z-scored within subject\n\n")
            elif task == "rrst":
                f.write("  Correct = 1 (upper), Incorrect = -1 (lower)\n")
                f.write("  stim_c = stimulus intensity, z-scored\n\n")

            # Summary table
            if r["summary"] is not None:
                f.write("Parameter Estimates:\n")
                f.write(r["summary"].to_string())
            else:
                f.write("(summary unavailable)\n")
            f.write("\n\n")

    print(f"  Report saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate HSSM report")
    parser.add_argument("--task", choices=["hrd", "rrst", "both"])
    parser.add_argument("--run-id", nargs="+",
                        help="Specific run ID(s)")
    parser.add_argument("--output", help="Output file path")
    args = parser.parse_args()

    # Gather results
    results = []
    if args.run_id:
        for rid in args.run_id:
            results.append(load_fitted(run_id=rid))
    elif args.task:
        tasks = ["hrd", "rrst"] if args.task == "both" else [args.task]
        for t in tasks:
            available = list_fitted(t)
            if available:
                results.append(load_fitted(run_id=available[0]))
            else:
                print(f"  No fitted models for {t.upper()}")
    else:
        # Default: most recent of both
        for t in ["hrd", "rrst"]:
            available = list_fitted(t)
            if available:
                results.append(load_fitted(run_id=available[0]))

    if not results:
        print("  No fitted models found.")
        return

    # Output path
    if args.output:
        output_path = Path(args.output)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = DOCS_DIR / f"hssm_report_{ts}.txt"

    print(f"\nGenerating report for {len(results)} model(s):")
    for r in results:
        print(f"  - {r['run_id']}")

    generate_report(results, output_path)
    print(f"\n====== Done! ======")


if __name__ == "__main__":
    main()
