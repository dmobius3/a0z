"""
jwst_speedup.py
===============

JWST early-galaxy free-fall collapse-time speedup applied to the Labbé
candidate sample (paper §3.4).

Setup (paper §3.4)
------------------

Under the framework's evolving acceleration scale a_0(z) = a_0(0) * E(z)
(eq. 3.1, anchored so E(0) = 1), the deep-MOND effective gravitational
acceleration at fixed Newtonian source

    g_eff = sqrt( g_N * a_0(z) )                                (§3.4)

scales as sqrt(E(z)). The corresponding free-fall collapse time

    t_ff propto 1 / sqrt( g_eff )

shortens as E(z)^(-1/4). Equivalently, the required star-formation
efficiency to assemble a given M_star at observed z is reduced from the
LambdaCDM value by

    epsilon_SF^MIT(z) = epsilon_SF^LCDM(z) / E(z)^(1/4).        (§3.4)

The framework speedup factor on collapse time is therefore E(z)^(1/4),
identical (by exponent identity) to the §3.2 v_flat enhancement and to the
§3.1 BTFR velocity ratio v(z)/v(0). To stay consistent with Modules 3-5,
this script evaluates the *anchored* ratio E(z)/E(0), so that
a_0(0) = a_0_SPARC = 1.20e-10 m/s^2 exactly. At z ~ 7-9 the un-normalized
and anchored E differ at the 5th significant figure; both round to the
same two-decimal-place display the paper quotes.

Sample
------

Six "central-value" massive candidates from Labbé et al. 2023
(Nature 616, 266) at log10(M_star / M_sun) > 10:

    id     z_phot
    11184  7.318
    38094  7.477
    2859   8.106
    13050  8.137
    14924  8.831
    35300  9.077

Per-candidate z_phot values are taken from the published Labbé candidate
table (the `sample_revision3_2207.12446.ecsv` file accompanying the paper)
at the catalog's full numerical precision; see
`../data/labbe_2023_candidates.csv` for the local copy.

Verification target (§3.4)
--------------------------

The paper's quantitative claim is the speedup-factor range:

    E(z)^(1/4) ∈ [1.92, 2.06]          (§3.4)

evaluated across the six z_phot values, with a corresponding
"approximately a factor of two" uniform reduction in epsilon_SF.

Imports
-------

    cosmology.PLANCK18  for E(z)
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import List, Sequence, Tuple

from cosmology import PLANCK18, Cosmology


# ---------------------------------------------------------------------------
# Anchored E(z) ratio (Modules 3-5 convention)
# ---------------------------------------------------------------------------

def E_eff(z, cosmo: Cosmology = PLANCK18) -> float:
    """E(z)/E(0): the anchored Friedmann ratio used throughout the paper's
    §3 onward, so that a_0(0) reduces to the SPARC-calibrated local value
    at z = 0 exactly. Mirrors btfr._E_eff."""
    return float(cosmo.E(z)) / float(cosmo.E(0.0))


def speedup_factor(z, cosmo: Cosmology = PLANCK18) -> float:
    """The §3.4 free-fall collapse-time speedup factor, E(z)^(1/4) under
    the anchored E(0) = 1 convention.

    By the §3.4 derivation:
        t_ff(z) / t_ff(0) = E(z)^(-1/4)
        speedup(z) = t_ff(0) / t_ff(z) = E(z)^(+1/4).
    """
    return E_eff(z, cosmo=cosmo) ** 0.25


def epsilon_reduction(z, cosmo: Cosmology = PLANCK18) -> float:
    """Multiplicative reduction of the required star-formation efficiency:
    epsilon_SF^MIT / epsilon_SF^LCDM = 1 / E(z)^(1/4)."""
    return 1.0 / speedup_factor(z, cosmo=cosmo)


# ---------------------------------------------------------------------------
# Labbé candidate loader
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LabbeCandidate:
    id: int
    z_phot: float
    log_mstar: float


# Default sample (full-precision z_phot values from the Labbé published
# candidate table `sample_revision3_2207.12446.ecsv`). The script accepts a
# `--data` argument (or environment override) pointing to a local copy.
DEFAULT_LABBE_SAMPLE: Tuple[LabbeCandidate, ...] = (
    LabbeCandidate(11184, 7.317981243133545, 10.181049329619398),
    LabbeCandidate(38094, 7.477300222506212, 10.886933911292022),
    LabbeCandidate( 2859, 8.105592642974775, 10.029126864461965),
    LabbeCandidate(13050, 8.136669296217242, 10.137980622560049),
    LabbeCandidate(14924, 8.830841454906686, 10.015011087250011),
    LabbeCandidate(35300, 9.077045913839543, 10.396785695245420),
)


def load_labbe_csv(path: str) -> List[LabbeCandidate]:
    """Parse `../data/labbe_2023_candidates.csv` per the schema in
    data/README.md."""
    out: List[LabbeCandidate] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            out.append(LabbeCandidate(
                id=int(row["id"]),
                z_phot=float(row["z_phot"]),
                log_mstar=float(row["log_mstar"]) if row.get("log_mstar") else float("nan"),
            ))
    return out


# ---------------------------------------------------------------------------
# Per-candidate table (data deposit row)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SpeedupRow:
    id: int
    z_phot: float
    E_eff: float                # E(z) / E(0) under the anchored convention
    speedup: float              # E_eff^(1/4)
    epsilon_factor: float       # 1 / speedup, the eps_SF reduction factor


def speedup_table(sample: Sequence[LabbeCandidate] = DEFAULT_LABBE_SAMPLE,
                  cosmo: Cosmology = PLANCK18) -> List[SpeedupRow]:
    rows: List[SpeedupRow] = []
    for c in sample:
        Eeff = E_eff(c.z_phot, cosmo=cosmo)
        s = Eeff ** 0.25
        rows.append(SpeedupRow(
            id=c.id,
            z_phot=c.z_phot,
            E_eff=Eeff,
            speedup=s,
            epsilon_factor=1.0 / s,
        ))
    return rows


# ---------------------------------------------------------------------------
# Verification driver
# ---------------------------------------------------------------------------

PAPER_RANGE_LOW = 1.92
PAPER_RANGE_HIGH = 2.06


def main():
    print("§3.4 JWST early-galaxy collapse-time speedup")
    print("=" * 60)
    print()
    print("Per-candidate speedup factors E(z)^(1/4) under anchored E(0)=1:")
    print(f"{'id':>6} {'z_phot':>10} {'E(z)':>10} {'E_eff':>10} "
          f"{'E^(1/4)':>10} {'1/E^(1/4)':>10}")
    rows = speedup_table()
    for r in rows:
        Ez = float(PLANCK18.E(r.z_phot))
        print(f"{r.id:6d} {r.z_phot:10.6f} {Ez:10.6f} {r.E_eff:10.6f} "
              f"{r.speedup:10.6f} {r.epsilon_factor:10.6f}")

    speedups = [r.speedup for r in rows]
    factors = [r.epsilon_factor for r in rows]
    s_min, s_max = min(speedups), max(speedups)
    eps_min, eps_max = min(factors), max(factors)

    print()
    print(f"Speedup range          : {s_min:.6f} -- {s_max:.6f}")
    print(f"  Rounded to 2 dp      : {round(s_min, 2)} -- {round(s_max, 2)}")
    print(f"  Paper §3.4 claim     : {PAPER_RANGE_LOW} -- {PAPER_RANGE_HIGH}")
    print()
    print(f"Epsilon-reduction range: {eps_min:.6f} -- {eps_max:.6f}")
    print(f"  Mean reduction       : {sum(factors)/len(factors):.6f}")
    print(f"  Paper §3.4 claim     : 'approximately a factor of two'")
    print()

    # ----- Pass / fail logic -----
    s_min_dp2 = round(s_min, 2)
    s_max_dp2 = round(s_max, 2)
    ok_low = s_min_dp2 == PAPER_RANGE_LOW
    ok_high = s_max_dp2 == PAPER_RANGE_HIGH

    # "approximately a factor of two" tolerance: each per-candidate epsilon
    # reduction should round to 0.5 (i.e. 0.45 <= 1/E^(1/4) <= 0.55) and the
    # mean reduction should round to 0.5 to 1 dp.
    ok_factor_two_each = all(round(f, 1) == 0.5 for f in factors)
    mean_factor = sum(factors) / len(factors)
    ok_factor_two_mean = round(mean_factor, 1) == 0.5

    print("Verification:")
    print(f"  E^(1/4) min @ 2 dp = {s_min_dp2}, paper = {PAPER_RANGE_LOW}: "
          f"{'PASS' if ok_low else 'FAIL'}")
    print(f"  E^(1/4) max @ 2 dp = {s_max_dp2}, paper = {PAPER_RANGE_HIGH}: "
          f"{'PASS' if ok_high else 'FAIL'}")
    print(f"  All per-candidate 1/E^(1/4) round to 0.5: "
          f"{'PASS' if ok_factor_two_each else 'FAIL'}")
    print(f"  Mean 1/E^(1/4) rounds to 0.5: "
          f"{'PASS' if ok_factor_two_mean else 'FAIL'}")

    print()
    print("=" * 60)
    if ok_low and ok_high and ok_factor_two_each and ok_factor_two_mean:
        print("ALL §3.4 NUMERICAL CLAIMS MATCH THE PAPER AT DISPLAYED PRECISION.")
    else:
        print("DISAGREEMENTS FOUND:")
        if not ok_low:
            print(f"  speedup min: paper={PAPER_RANGE_LOW}, computed={s_min_dp2}")
        if not ok_high:
            print(f"  speedup max: paper={PAPER_RANGE_HIGH}, computed={s_max_dp2}")
        if not ok_factor_two_each:
            print(f"  per-candidate factor-of-two: not all round to 0.5")
        if not ok_factor_two_mean:
            print(f"  mean factor-of-two: rounds to {round(mean_factor,1)}, not 0.5")


if __name__ == "__main__":
    main()
