"""Single-command pass/fail verification harness for the a0z pipeline.

Runs every analysis module's `main()` self-check, plus `make_tables.py`'s
cell-by-cell verification, and aggregates the results.

Usage
-----
    python verify.py             # run every module, aggregate pass/fail
    python verify.py --quiet     # suppress per-module output, summary only
    python verify.py --module framework  # run a single module's self-check

Exit code 0 = every module green.
Exit code 1 = any module failed (per-module non-zero exit, "FAIL" in stdout,
              "match: False" patterns, or "Traceback" in stderr).

The harness shells out to each script rather than importing API points,
so a future change to a module's docstring or main() output format
cannot quietly invalidate the gate.
"""

import argparse
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from typing import List, Optional


SEED = 42  # locked for any forthcoming Monte Carlo (Übler etc.)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Modules in the order they should be run. Tier 1 first (no dependencies),
# then prediction modules, then make_tables.py last (depends on the rest).
MODULES = [
    "cosmology",
    "framework",
    "combinatorial_baseline",
    "btfr",
    "rotation_curves",
    "lensing",
    "jwst_speedup",
    "ubler_sigma_tension",
    "cmb_leakage",
    "make_tables",
]

# Output-string patterns that indicate failure even when a module's
# main() exits 0. Any module is considered failed if any of these
# substrings appears in its stdout (case-insensitive for some).
FAILURE_PATTERNS = [
    re.compile(r"\bFAIL\b"),
    re.compile(r"Match\s*\(\d+\s*dp\)\s*:\s*False"),
    re.compile(r"match:\s*False", re.IGNORECASE),
    re.compile(r"Error\b", re.IGNORECASE),
]


@dataclass
class ModuleResult:
    name: str
    exit_code: int
    stdout: str
    stderr: str
    ok: bool
    failure_reason: Optional[str]


def _summary_line(stdout: str) -> str:
    """Pull a 1-line summary from a module's stdout for the harness report.

    Looks for the canonical 'ALL ... MATCH' / 'All ... match: True'
    style closing line, falling back to the last non-empty line.
    """
    candidates = [
        re.compile(r"^ALL .* MATCH .*$", re.MULTILINE | re.IGNORECASE),
        re.compile(r"^All .*claims match.*$", re.MULTILINE | re.IGNORECASE),
        re.compile(r"^ALL CELLS IN ALL TABLES MATCH .*$",
                   re.MULTILINE | re.IGNORECASE),
        re.compile(r"^Match\s*\(\d+\s*dp\)\s*:\s*True\b.*$", re.MULTILINE),
    ]
    for pat in candidates:
        m = pat.search(stdout)
        if m:
            return m.group(0).strip()
    # Fallback: last non-empty line
    for line in reversed(stdout.splitlines()):
        if line.strip():
            return line.strip()
    return ""


def run_module(name: str) -> ModuleResult:
    """Run scripts/<name>.py via subprocess and classify pass/fail."""
    path = os.path.join(SCRIPT_DIR, f"{name}.py")
    proc = subprocess.run(
        [sys.executable, path],
        capture_output=True, text=True,
        cwd=SCRIPT_DIR,
    )
    failure_reason: Optional[str] = None

    if "Traceback" in proc.stderr or "Error" in proc.stderr:
        failure_reason = "stderr contains traceback or error"
    elif proc.returncode != 0:
        failure_reason = f"non-zero exit code {proc.returncode}"
    else:
        for pat in FAILURE_PATTERNS:
            if pat.search(proc.stdout):
                failure_reason = f"stdout matches failure pattern {pat.pattern!r}"
                break

    return ModuleResult(
        name=name,
        exit_code=proc.returncode,
        stdout=proc.stdout,
        stderr=proc.stderr,
        ok=(failure_reason is None),
        failure_reason=failure_reason,
    )


def run(modules: List[str], quiet: bool) -> int:
    print()
    print("=" * 72)
    print(f" a0z verification harness    seed={SEED}")
    print("=" * 72)
    print()

    results: List[ModuleResult] = []
    for name in modules:
        r = run_module(name)
        results.append(r)
        flag = "OK  " if r.ok else "FAIL"
        summary = _summary_line(r.stdout) if r.ok else (r.failure_reason or "")
        print(f"  [{flag}]  {name:<26} {summary}")

    n_pass = sum(1 for r in results if r.ok)
    n_total = len(results)
    print()
    print("-" * 72)
    print(f" {n_pass}/{n_total} modules pass")
    print("-" * 72)
    print()

    failures = [r for r in results if not r.ok]
    if failures:
        if not quiet:
            for r in failures:
                print(f"=== FAILED: {r.name} ===")
                print(f"  reason   : {r.failure_reason}")
                print(f"  exit code: {r.exit_code}")
                if r.stderr.strip():
                    print(f"  stderr   :")
                    for line in r.stderr.splitlines():
                        print(f"    {line}")
                if r.stdout.strip():
                    # Only print the tail of stdout to keep noise down.
                    tail = r.stdout.splitlines()[-15:]
                    print(f"  stdout (tail):")
                    for line in tail:
                        print(f"    {line}")
                print()
        return 1

    print("All modules report green against the lean PRD draft of the paper.")
    print()
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n", 1)[0],
    )
    parser.add_argument(
        "--module",
        action="append",
        default=None,
        help="Run only the named module(s). Repeat for multiple.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress per-failure stdout/stderr dumps; show summary only.",
    )
    args = parser.parse_args()

    if args.module:
        unknown = [m for m in args.module if m not in MODULES]
        if unknown:
            print(f"unknown module(s): {unknown}", file=sys.stderr)
            print(f"available: {MODULES}", file=sys.stderr)
            return 2
        modules = args.module
    else:
        modules = MODULES

    return run(modules, quiet=args.quiet)


if __name__ == "__main__":
    sys.exit(main())
