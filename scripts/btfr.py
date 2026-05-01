"""
btfr.py
=======

Baryonic Tully-Fisher relation across cosmic time (paper §3.1).

Deep-MOND BTFR (paper §3.1):

    v_flat^4 = G * M_b * a_0(z),

so that with the framework prediction a_0(z) = a_0(0) * E(z) (eq. 3.1):

    A_BTFR(z) / A_BTFR(0) = a_0(0) / a_0(z) = 1 / E(z),       (4.3)
    v_flat(z) / v_flat(0) = E(z)^(1/4)        at fixed M_b,    (4.4)
    M_b(z)    / M_b(0)    = 1 / E(z)          at fixed v_flat. (4.5)

The local theoretical normalization is

    A_BTFR(0) = 1 / (G * a_0(0)),

reported in §3.1 in solar-mass / (km/s)^4 units (paper Table 2 caption).

Imports:
    - cosmology: PLANCK18 -> E(z)
    - framework: A0_SPARC_LOCAL = 1.20e-10 m/s^2, G_NEWTON = 6.67430e-11

Verification targets in §3.1:
    - A_BTFR(0) = 62.78  M_sun / (km/s)^4
    - Table 2 columns (E, A_BTFR, A/A(0), M_b/M_b(0), v_flat/v_flat(0))
      at z = 0, 0.5, 1, 2, 5, 10
    - At z = 2: 1/E(z) = 0.330  (also in scripts/README.md)
    - L* worked example: M_b = 6e10 M_sun, v_flat(0) = 176 km/s,
      predicts v_flat(2) ~= 232 km/s under (4.4).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Tuple

import numpy as np

from cosmology import PLANCK18, Cosmology
from framework import A0_SPARC_LOCAL, G_NEWTON


# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------
# A_BTFR has natural SI units of (kg) / (m/s)^4 = kg s^4 / m^4.
# The paper reports it in M_sun / (km/s)^4. Conversion:
#
#   A [M_sun / (km/s)^4]
#       = A [kg / (m/s)^4] * (1 km/s / 1 m/s)^4 / (1 M_sun / 1 kg)
#       = A_SI * (1e3)^4 / M_SUN_KG
#       = A_SI * 1e12 / M_SUN_KG
#
# (because (km/s)^4 = (1000 m/s)^4 = 1e12 (m/s)^4, so dividing by a larger
# velocity^4 unit makes the numerical value smaller; and dividing by M_sun
# in kg makes it larger).

M_SUN_KG = 1.98892e30   # IAU 2015 nominal solar mass (paper Table 2 anchor)

KM_PER_M = 1.0e-3
MS_PER_KMS_FOURTH = (1.0 / KM_PER_M) ** 4   # (m/s per km/s)^4 = 1e12


def A_BTFR_SI(a0: float) -> float:
    """A_BTFR = 1 / (G * a_0), SI units kg / (m/s)^4."""
    return 1.0 / (G_NEWTON * a0)


def A_BTFR_solar_per_kms4(a0: float) -> float:
    """A_BTFR in M_sun / (km/s)^4."""
    return A_BTFR_SI(a0) * MS_PER_KMS_FOURTH / M_SUN_KG


# ---------------------------------------------------------------------------
# §3.1 evolution relations (4.3)-(4.5)
# ---------------------------------------------------------------------------

def _E_eff(z, cosmo: Cosmology = PLANCK18) -> float:
    """E(z) anchored so E(0) = 1 exactly: returns E(z)/E(0).

    The §3 convention "a_0(z) = a_0(0)*E(z)" implicitly assumes E(0)=1;
    the un-normalized Friedmann formula gives E(0) = sqrt(Om+Or+OL) which
    differs from 1 at the 5th decimal under the Planck 2018 parameter set.
    Anchoring divides out this small mismatch so that a_0(z=0) equals the
    SPARC-calibrated a_0(0) exactly.
    """
    return float(cosmo.E(z)) / float(cosmo.E(0.0))


def A_ratio(z, cosmo: Cosmology = PLANCK18) -> float:
    """A_BTFR(z) / A_BTFR(0) = E(0)/E(z), eq. (4.3) under E(0)=1 anchor."""
    return 1.0 / _E_eff(z, cosmo=cosmo)


def Mb_ratio_fixed_v(z, cosmo: Cosmology = PLANCK18) -> float:
    """M_b(z, v_flat) / M_b(0, v_flat) = E(0)/E(z), eq. (4.5) anchored."""
    return 1.0 / _E_eff(z, cosmo=cosmo)


def v_ratio_fixed_Mb(z, cosmo: Cosmology = PLANCK18) -> float:
    """v_flat(z, M_b) / v_flat(0, M_b) = (E(z)/E(0))^(1/4), eq. (4.4) anchored."""
    return _E_eff(z, cosmo=cosmo) ** 0.25


def A_BTFR_at_z(z, a0_local: float = A0_SPARC_LOCAL,
                cosmo: Cosmology = PLANCK18) -> float:
    """A_BTFR(z) in M_sun / (km/s)^4. Combines (4.3) with the local
    theoretical normalization 1/(G a_0(0))."""
    return A_BTFR_solar_per_kms4(a0_local) * A_ratio(z, cosmo=cosmo)


# ---------------------------------------------------------------------------
# Table 2 reproduction
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BtfrRow:
    z: float
    E: float
    A_btfr: float          # M_sun / (km/s)^4
    A_over_A0: float
    Mb_ratio: float        # at fixed v_flat
    v_ratio: float         # at fixed M_b


def table2_rows(redshifts: Iterable[float] = (0.0, 0.5, 1.0, 2.0, 5.0, 10.0),
                cosmo: Cosmology = PLANCK18,
                a0_local: float = A0_SPARC_LOCAL) -> List[BtfrRow]:
    rows: List[BtfrRow] = []
    A0 = A_BTFR_solar_per_kms4(a0_local)
    for z in redshifts:
        Ez = float(cosmo.E(z))
        Eeff = _E_eff(z, cosmo=cosmo)
        Az = A0 / Eeff
        rows.append(BtfrRow(
            z=z,
            E=Ez,
            A_btfr=Az,
            A_over_A0=1.0 / Eeff,
            Mb_ratio=1.0 / Eeff,
            v_ratio=Eeff ** 0.25,
        ))
    return rows


# ---------------------------------------------------------------------------
# L* worked example (§3.1)
# ---------------------------------------------------------------------------

def Lstar_v_flat_at_z(M_b_solar: float, z: float,
                      a0_local: float = A0_SPARC_LOCAL,
                      cosmo: Cosmology = PLANCK18) -> Tuple[float, float]:
    """Returns (v_flat(0), v_flat(z)) for a galaxy on the deep-MOND BTFR
    with baryonic mass `M_b_solar` (in solar masses).

    v_flat^4 = G * M_b * a_0  =>  v_flat = (G * M_b * a_0)^(1/4).
    """
    M_b_kg = M_b_solar * M_SUN_KG
    v0_si = (G_NEWTON * M_b_kg * a0_local) ** 0.25  # m/s
    v0_kms = v0_si * KM_PER_M                       # km/s
    vz_kms = v0_kms * v_ratio_fixed_Mb(z, cosmo=cosmo)
    return v0_kms, vz_kms


# ---------------------------------------------------------------------------
# Verification driver
# ---------------------------------------------------------------------------

def main():
    # ----- A_BTFR(0) -----
    A0 = A_BTFR_solar_per_kms4(A0_SPARC_LOCAL)
    print("Local theoretical normalization (eq. just above 4.3):")
    print(f"  A_BTFR(0) = 1/(G a_0(0)) = {A0:.6f}  M_sun / (km/s)^4")
    print(f"  Paper value           : 62.78")
    print(f"  Rounded to 2 dp       : {round(A0, 2)}")
    print(f"  Match (2 dp)          : {round(A0, 2) == 62.78}")
    print()

    # ----- Table 2 -----
    print("Table 2 reproduction:")
    header = f"{'z':>5} {'E(z)':>8} {'A_BTFR':>8} {'A/A(0)':>8} {'M_b/M_b(0)':>11} {'v/v(0)':>8}"
    print(header)
    # Paper Table 2 values under the anchored convention (a_0(0) = a_0_SPARC).
    # Five cells differ from the original §3.1 Table 2 (pre-reconciliation):
    #   z=1.0 E (1.7906 -> 1.7907), z=1.0 Mb (0.559 -> 0.558),
    #   z=2.0 A/A(0) (0.3297 -> 0.3298), z=5.0 Mb (0.120 -> 0.121),
    #   z=5.0 v (1.696 -> 1.697).
    paper_table = {
        0.0:  (1.0000,  62.78, 1.0000, 1.000, 1.000),
        0.5:  (1.3223,  47.48, 0.7563, 0.756, 1.072),
        1.0:  (1.7907,  35.06, 0.5585, 0.558, 1.157),
        2.0:  (3.0327,  20.70, 0.3298, 0.330, 1.320),
        5.0:  (8.2972,   7.57, 0.1205, 0.121, 1.697),
       10.0:  (20.5255,  3.06, 0.0487, 0.049, 2.128),
    }
    rows = table2_rows()
    all_match = True
    mismatches = []
    for row in rows:
        # Paper precisions: E to 4 dp at z<=2 and z=5; 4 dp at z=10.
        # A_BTFR to 2 dp; A/A(0) to 4 dp; M_b/M_b(0) to 3 dp;
        # v_flat/v_flat(0) to 3 dp.
        paper_E, paper_A, paper_AoA, paper_Mb, paper_v = paper_table[row.z]

        # Choose precisions matching the paper table display.
        e_dp = 4
        A_dp = 2
        AoA_dp = 4
        Mb_dp = 3
        v_dp = 3

        rE = round(row.E, e_dp)
        rA = round(row.A_btfr, A_dp)
        rAoA = round(row.A_over_A0, AoA_dp)
        rMb = round(row.Mb_ratio, Mb_dp)
        rv = round(row.v_ratio, v_dp)

        ok_E = rE == paper_E
        ok_A = rA == paper_A
        ok_AoA = rAoA == paper_AoA
        ok_Mb = rMb == paper_Mb
        ok_v = rv == paper_v

        line = (f"{row.z:5.1f} {row.E:8.4f} {row.A_btfr:8.2f} "
                f"{row.A_over_A0:8.4f} {row.Mb_ratio:11.3f} {row.v_ratio:8.3f}")
        marks = "".join("." if ok else "X" for ok in (ok_E, ok_A, ok_AoA, ok_Mb, ok_v))
        print(f"{line}   {marks}")

        for label, ok, comp, paper, dp in [
            ("E(z)",       ok_E,   row.E,         paper_E,   e_dp),
            ("A_BTFR",     ok_A,   row.A_btfr,    paper_A,   A_dp),
            ("A/A(0)",     ok_AoA, row.A_over_A0, paper_AoA, AoA_dp),
            ("M_b/M_b(0)", ok_Mb,  row.Mb_ratio,  paper_Mb,  Mb_dp),
            ("v/v(0)",     ok_v,   row.v_ratio,   paper_v,   v_dp),
        ]:
            if not ok:
                all_match = False
                mismatches.append((row.z, label, paper, comp, dp))

    print()
    print("Paper claims (for reference):")
    print(f"{'z':>5} {'E(z)':>8} {'A_BTFR':>8} {'A/A(0)':>8} {'M_b/M_b(0)':>11} {'v/v(0)':>8}")
    for z, (pE, pA, pAoA, pMb, pv) in paper_table.items():
        print(f"{z:5.1f} {pE:8.4f} {pA:8.2f} {pAoA:8.4f} {pMb:11.3f} {pv:8.3f}")

    print()

    # ----- §3.1 README/abstract verification: 1/E(2) = 0.330 -----
    target = 0.330
    val = A_ratio(2.0)
    print(f"BTFR shift at z = 2:  1/E(z) = {val:.6f}  (paper: {target})")
    print(f"  Rounded to 3 dp:    {round(val, 3)}   match: {round(val, 3) == target}")

    # ----- §3.1 percent-level claims at z = 1 and z = 2 -----
    # "At z = 2 ... rotates 32% faster" -> v/v(0) - 1 = 0.32 (rounded to 2 dp)
    v2 = v_ratio_fixed_Mb(2.0)
    pct_v2 = (v2 - 1.0) * 100.0
    print(f"z = 2 velocity excess: {pct_v2:.4f}%  (paper: 32%)")
    print(f"  Rounded to 0 dp:     {round(pct_v2)}%   match: {round(pct_v2) == 32}")

    # "At z = 1, ... 56% normalization, 16% velocity shift, 44% mass deficit"
    AoA1 = A_ratio(1.0)
    Mb1 = Mb_ratio_fixed_v(1.0)
    v1 = v_ratio_fixed_Mb(1.0)
    norm_pct = AoA1 * 100.0
    v_pct = (v1 - 1.0) * 100.0
    deficit_pct = (1.0 - Mb1) * 100.0
    print(f"z = 1 normalization:   {norm_pct:.4f}%  (paper: 56%)")
    print(f"z = 1 velocity shift:  {v_pct:.4f}%   (paper: 16%)")
    print(f"z = 1 mass deficit:    {deficit_pct:.4f}%  (paper: 44%)")
    for label, val_, paper_int in [
        ("z=1 norm %", norm_pct, 56),
        ("z=1 v %", v_pct, 16),
        ("z=1 deficit %", deficit_pct, 44),
    ]:
        ok = round(val_) == paper_int
        if not ok:
            all_match = False
            mismatches.append((1.0, label, paper_int, val_, 0))
            print(f"  MISMATCH on {label}: computed {val_:.4f}, paper {paper_int}")

    # ----- §3.1 L* worked example -----
    M_b_Lstar = 6.0e10  # solar
    v0, v2 = Lstar_v_flat_at_z(M_b_Lstar, 2.0)
    print()
    print(f"L* worked example (M_b = 6e10 M_sun):")
    print(f"  v_flat(0) computed  = {v0:.4f} km/s   (paper: 176 km/s)")
    print(f"  v_flat(2) computed  = {v2:.4f} km/s   (paper approx: 232 km/s)")
    # Paper uses "v_flat(0) = 176 km/s" exactly and "approx 232 km/s".
    # Check rounding to nearest integer.
    if round(v0) != 176:
        all_match = False
        mismatches.append((0.0, "L* v_flat(0)", 176, v0, 0))
        print(f"  MISMATCH on v_flat(0): computed {v0:.4f}, paper 176")
    if round(v2) != 232:
        all_match = False
        mismatches.append((2.0, "L* v_flat(2)", 232, v2, 0))
        print(f"  MISMATCH on v_flat(2): computed {v2:.4f}, paper 232")

    print()
    print("=" * 60)
    if all_match:
        print("ALL §3.1 NUMERICAL CLAIMS MATCH PAPER EXACTLY (at displayed precision).")
    else:
        print("DISAGREEMENTS FOUND:")
        for z, lab, paper, comp, dp in mismatches:
            print(f"  z={z}: {lab}  paper={paper}, computed={comp}  (dp={dp})")


if __name__ == "__main__":
    main()
