"""
combinatorial_baseline.py
=========================

Module 9 of the a0z analysis pipeline. Computes the combinatorial sparsity
baseline against which the §2.5 well-pair match (a_0/(c H_0) recovered from
the (13/120, 34/120) Fibonacci-well pair) is calibrated.

Paper claim (§2.5)
------------------
"Of the 7,021 unordered nonzero phase-position pairs (k_1, k_2)/120 with
k_1, k_2 in {1, ..., 119}, only 24 produce a ratio C(k_1/120)/C(k_2/120)
within 1% of the observed 0.1832 [post-reconciliation: 0.1833], a fraction
of 0.342%. Modulo the reflection symmetry C(k) = C(120-k), these collapse
to 6 unique phase-operator value pairs. Within the framework's selected
Fibonacci-well subset {13, 21, 34, 55, 60}/120, exactly one of the ten
possible unordered well-pairs matches: the framework's (13, 34). The match
is therefore sparse on the full 120-domain (one in ~300 pairs) and unique
within the topologically-selected Fibonacci-well structure (one in 10)."

Conventions used here
---------------------
- Phase operator: C(Theta) = 2 sin^2(pi Theta), imported from framework.py.
- Domain: k in {1, ..., 119} on the 120-domain (k = 0 excluded since C(0) = 0
  would make the ratio undefined, matching the paper's "nonzero" wording).
- Unordered pair: each (k_1, k_2) with k_1 < k_2 counted once. The number of
  such pairs is C(119, 2) = 7,021.
- "Ratio within 1%": each unordered pair (k_1, k_2) yields two ratios,
  r and 1/r; only one of them lies near the small target 0.1833. We take the
  smaller-of-two (i.e. the one < 1) and test |r - 0.1833| / 0.1833 <= 0.01.
- "Modulo reflection symmetry": C(k/120) = C((120-k)/120), so the canonical
  representative of k is min(k, 120-k). Two unordered pairs collapse if their
  canonical-representative sets coincide.
- Fibonacci-well subset: {13, 21, 34, 55, 60}/120, giving C(5, 2) = 10
  unordered pairs.
- Observed Milgrom ratio: 0.1833 (post §2 reconciliation; the paper text
  printed 0.1832 prior to the reconciliation pass).

Verification target (per scripts/README.md)
-------------------------------------------
"Combinatorial sparsity: 24 of 7,021 phase-position pairs match within 1%;
1 of 10 within Fibonacci-well subset"

Run
---
    python combinatorial_baseline.py
        -> prints a verification table and exits 0 on full agreement.
"""

from __future__ import annotations

import argparse
import itertools
import sys
from dataclasses import dataclass
from typing import List, Tuple

from framework import C_phase, DOMAIN_SIZE


# ---------------------------------------------------------------------------
# Paper inputs (post §2 reconciliation)
# ---------------------------------------------------------------------------

OBSERVED_RATIO = 0.1833           # a_0_obs / (c H_0_obs), post-reconciliation
TOLERANCE_FRAC = 0.01             # "within 1%"

FIBONACCI_WELLS = (13, 21, 34, 55, 60)
FRAMEWORK_PAIR = (13, 34)         # the well pair matched by the framework


# ---------------------------------------------------------------------------
# Core combinatorial computation
# ---------------------------------------------------------------------------

def smaller_ratio(k1: int, k2: int, N: int = DOMAIN_SIZE) -> float:
    """
    For an unordered pair (k1, k2), return the smaller of C(k1/N)/C(k2/N) and
    C(k2/N)/C(k1/N). This is the ratio in (0, 1]; it is the one that can lie
    near 0.1833. The other ratio is its reciprocal, near ~5.45, which never
    matches the small target by construction.
    """
    c1 = float(C_phase(k1 / N))
    c2 = float(C_phase(k2 / N))
    r = c1 / c2 if c1 < c2 else c2 / c1
    return r


def matches_within(
    k1: int,
    k2: int,
    target: float = OBSERVED_RATIO,
    tol: float = TOLERANCE_FRAC,
    N: int = DOMAIN_SIZE,
) -> bool:
    """True iff the smaller-of-two ratio for (k1, k2) lies within `tol` of `target`."""
    r = smaller_ratio(k1, k2, N=N)
    return abs(r - target) / target <= tol


def canonical_under_reflection(
    k1: int, k2: int, N: int = DOMAIN_SIZE
) -> Tuple[int, int]:
    """
    Reduce each index k to min(k, N-k), then sort. Two pairs that map to the
    same canonical tuple share the same (C-value, C-value) set, since
    C(k/N) = C((N-k)/N).
    """
    a = min(k1, N - k1)
    b = min(k2, N - k2)
    return tuple(sorted((a, b)))


@dataclass(frozen=True)
class CombinatorialResult:
    total_pairs: int                       # 7021
    matching_pairs: int                    # 24
    matching_fraction_percent: float       # 0.342...
    unique_under_reflection: int           # 6
    fibonacci_total_pairs: int             # 10
    fibonacci_matching_pairs: int          # 1
    fibonacci_match_list: List[Tuple[int, int]]
    sparse_one_in_n: float                 # ~292.5 (i.e. ~300)


def run_combinatorial(
    N: int = DOMAIN_SIZE,
    target: float = OBSERVED_RATIO,
    tol: float = TOLERANCE_FRAC,
    fib_wells: Tuple[int, ...] = FIBONACCI_WELLS,
) -> CombinatorialResult:
    """
    Enumerate all unordered nonzero pairs on the N-domain, count the ones
    whose smaller-of-two phase-operator ratio lies within `tol` of `target`,
    and report the §2.5 statistics.
    """

    indices = list(range(1, N))                 # nonzero: exclude k = 0
    all_pairs = list(itertools.combinations(indices, 2))

    matching = [(k1, k2) for k1, k2 in all_pairs if matches_within(k1, k2, target, tol, N)]

    canonical_set = {canonical_under_reflection(k1, k2, N) for k1, k2 in matching}

    fib_pairs = list(itertools.combinations(fib_wells, 2))
    fib_match = [(k1, k2) for k1, k2 in fib_pairs if matches_within(k1, k2, target, tol, N)]

    return CombinatorialResult(
        total_pairs=len(all_pairs),
        matching_pairs=len(matching),
        matching_fraction_percent=100.0 * len(matching) / len(all_pairs),
        unique_under_reflection=len(canonical_set),
        fibonacci_total_pairs=len(fib_pairs),
        fibonacci_matching_pairs=len(fib_match),
        fibonacci_match_list=fib_match,
        sparse_one_in_n=len(all_pairs) / len(matching) if matching else float("inf"),
    )


# ---------------------------------------------------------------------------
# Verification against the §2.5 paper claims
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CheckRow:
    label: str
    paper: str
    computed: str
    matches: bool


def run_verification() -> Tuple[List[CheckRow], bool]:
    """Reproduce every §2.5 numerical claim and check it against the paper."""

    res = run_combinatorial()

    # 1. Total unordered nonzero pairs on the 120-domain.
    total_match = res.total_pairs == 7021

    # 2. Pairs matching within 1% of 0.1833.
    n_match = res.matching_pairs == 24

    # 3. Fraction (paper rounds 0.34183...% to 0.342%).
    fraction_str_3dp = f"{res.matching_fraction_percent:.3f}"
    fraction_match = fraction_str_3dp == "0.342"

    # 4. Unique phase-operator value pairs modulo reflection k -> 120-k.
    refl_match = res.unique_under_reflection == 6

    # 5. Fibonacci subset has 10 unordered pairs.
    fib_total_match = res.fibonacci_total_pairs == 10

    # 6. Exactly one of those ten pairs matches.
    fib_one_match = res.fibonacci_matching_pairs == 1

    # 7. The single Fibonacci match is the framework's (13, 34).
    fib_is_framework = (
        res.fibonacci_matching_pairs == 1
        and tuple(sorted(res.fibonacci_match_list[0])) == FRAMEWORK_PAIR
    )

    # 8. "Sparse on the full 120-domain (one in ~300 pairs)".
    # 7021 / 24 = 292.5..., paper rounds to ~300.
    sparse_match = round(res.sparse_one_in_n, -2) == 300.0

    rows = [
        CheckRow(
            "Total unordered nonzero pairs (k1, k2)/120",
            "7,021",
            f"{res.total_pairs}",
            total_match,
        ),
        CheckRow(
            "Pairs with smaller ratio within 1% of 0.1833",
            "24",
            f"{res.matching_pairs}",
            n_match,
        ),
        CheckRow(
            "Matching fraction",
            "0.342%",
            f"{fraction_str_3dp}%",
            fraction_match,
        ),
        CheckRow(
            "Unique pairs modulo reflection C(k) = C(120-k)",
            "6",
            f"{res.unique_under_reflection}",
            refl_match,
        ),
        CheckRow(
            "Fibonacci-well subset unordered pairs",
            "10",
            f"{res.fibonacci_total_pairs}",
            fib_total_match,
        ),
        CheckRow(
            "Fibonacci-well pairs matching within 1%",
            "1",
            f"{res.fibonacci_matching_pairs}",
            fib_one_match,
        ),
        CheckRow(
            "Identity of the Fibonacci-well match",
            "(13, 34)",
            f"{tuple(sorted(res.fibonacci_match_list[0])) if res.fibonacci_match_list else None}",
            fib_is_framework,
        ),
        CheckRow(
            "Sparsity (1 in N pairs)",
            "~300",
            f"~{round(res.sparse_one_in_n):d}",
            sparse_match,
        ),
    ]

    all_pass = all(r.matches for r in rows)
    return rows, all_pass


def print_table(rows: List[CheckRow]) -> None:
    label_w = max(len(r.label) for r in rows)
    paper_w = max(len(r.paper) for r in rows)
    comp_w = max(len(r.computed) for r in rows)
    label_w = max(label_w, len("Quantity"))
    paper_w = max(paper_w, len("Paper §2.5"))
    comp_w = max(comp_w, len("Computed"))

    print(f"{'Quantity':<{label_w}}  {'Paper §2.5':<{paper_w}}  {'Computed':<{comp_w}}  Match")
    print("-" * (label_w + paper_w + comp_w + 14))
    for r in rows:
        ok = "OK" if r.matches else "FAIL"
        print(f"{r.label:<{label_w}}  {r.paper:<{paper_w}}  {r.computed:<{comp_w}}  {ok}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument(
        "--target",
        type=float,
        default=OBSERVED_RATIO,
        help="Observed Milgrom ratio (default 0.1833 post-§2 reconciliation).",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=TOLERANCE_FRAC,
        help="Fractional tolerance for 'within 1%%' (default 0.01).",
    )
    parser.add_argument(
        "--list-matches",
        action="store_true",
        help="Print all 24 matching unordered pairs with their ratios.",
    )
    args = parser.parse_args()

    rows, all_pass = run_verification()
    print_table(rows)

    if args.list_matches:
        print()
        print("All matching unordered pairs (smaller-of-two ratio within 1% of target):")
        target = args.target
        tol = args.tol
        pairs = []
        for k1, k2 in itertools.combinations(range(1, DOMAIN_SIZE), 2):
            r = smaller_ratio(k1, k2)
            if abs(r - target) / target <= tol:
                pairs.append((k1, k2, r, canonical_under_reflection(k1, k2)))
        pairs.sort(key=lambda t: (t[3], t[0], t[1]))
        for k1, k2, r, can in pairs:
            print(f"  ({k1:3d}, {k2:3d})   smaller ratio = {r:.6f}   canonical = {can}")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
