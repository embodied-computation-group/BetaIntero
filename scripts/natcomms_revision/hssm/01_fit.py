#!/usr/bin/env python
"""
Fit HSSM DDM models and save results to disk.

Usage:
    python 01_fit.py --mode debug --task hrd
    python 01_fit.py --mode full --task both --parallel
    python 01_fit.py --mode quick --task rrst
"""

import argparse
import subprocess
import sys
from pathlib import Path

from hssm_utils import (
    MODE_SETTINGS, TASK_LOADERS, TASK_FORMULAS,
    build_model, fit_model, save_fitted,
)


def main():
    parser = argparse.ArgumentParser(description="Fit HSSM DDM models")
    parser.add_argument(
        "--mode", choices=["debug", "quick", "full"], default="debug",
        help="debug=1chain x 100, quick=4x500, full=4x1000",
    )
    parser.add_argument(
        "--task", choices=["hrd", "rrst", "both"], default="both",
    )
    parser.add_argument(
        "--parallel", action="store_true",
        help="Fit HRD and RRST in parallel subprocesses",
    )
    args = parser.parse_args()

    # Parallel dispatch
    if args.parallel and args.task == "both":
        script = str(Path(__file__).resolve())
        procs = []
        for task in ("hrd", "rrst"):
            cmd = [sys.executable, script, "--mode", args.mode, "--task", task]
            print(f"  Spawning: {' '.join(cmd)}")
            procs.append(subprocess.Popen(
                cmd, stdout=sys.stdout, stderr=sys.stderr,
            ))
        codes = [p.wait() for p in procs]
        if any(c != 0 for c in codes):
            print(f"  WARNING: child exits {codes}")
        sys.exit(max(codes))

    # Sequential
    tasks = ["hrd", "rrst"] if args.task == "both" else [args.task]
    cfg = MODE_SETTINGS[args.mode]

    for task in tasks:
        print(f"\n{'=' * 60}")
        print(f"Fitting {task.upper()} Model")
        print(f"Mode: {args.mode}  |  Sampler: {cfg['sampler']}")
        print(f"{cfg['chains']} chains x {cfg['draws']} draws "
              f"(tune={cfg['tune']})")
        print("=" * 60)

        formulas = TASK_FORMULAS[task]
        for param, formula in formulas.items():
            print(f"  {param}: {formula}")
        print()

        data = TASK_LOADERS[task]()
        model = build_model(task, data)
        fit_model(model, cfg)

        # Print summary
        try:
            summary = model.summary()
            print(f"\n{task.upper()} Fixed Effects:")
            fixed = summary[~summary.index.str.contains("participant_id")]
            print(fixed.to_string())
        except Exception:
            pass

        run_id = save_fitted(model, data, task, args.mode, cfg)
        print(f"  Run ID: {run_id}")

    print(f"\n====== Done! ======")


if __name__ == "__main__":
    main()
