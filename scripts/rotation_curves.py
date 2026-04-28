"""
rotation_curves.py
==================

Rotation-curve predictions across cosmic time (paper §5).

Key relations:

  - MOND transition radius (eq. 5.1):
        r_M(z) = sqrt( G * M_b / a_0(z) )
  - Fractional shifts (eqs. 5.2, 5.3):
        r_M(z) / r_M(0)        = 1 / sqrt(E(z))   (mass-independent)
        v_flat(z) / v_flat(0)  = E(z)^(1/4)       (mass-independent, deep MOND)

The asymptotic flat velocity is the BTFR statement (4.1):
        v_flat = ( G * M_b * a_0(z) )^(1/4).

Convention (carried over from Module 3 / btfr.py):
  a_0(z) is anchored so that a_0(z=0) equals the SPARC value exactly:
        a_0(z) = A0_SPARC_LOCAL * cosmo.E(z) / cosmo.E(0).
  The §5 formulas inherit this anchoring through their a_0(z) factor;
  in particular, r_M and v_flat at z=0 reduce to the SPARC-calibrated
  values and the fractional shifts reduce to powers of E(z)/E(0).

Imports:
  - cosmology: PLANCK18 (Planck 2018 flat LCDM)
  - framework: A0_SPARC_LOCAL = 1.20e-10 m/s^2, G_NEWTON = 6.67430e-11
  - btfr:      M_SUN_KG, KM_PER_M, _E_eff (anchored E(z)/E(0))

§5 paper claims this script verifies:
  - §5.1 fractional r_M shifts: 0.574 (z=2), 0.347 (z=5), 0.221 (z=10)
  - §5.2 Table 3 archetypes (Dwarf / Sub-L* / L* / Giant) at z=0 and z=2
  - §5.3 L* asymptotes 176 km/s (z=0), 232 km/s (z=2)
  - Caption shifts: r_M(2)/r_M(0) = 0.574, v_flat(2)/v_flat(0) = 1.320
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Tuple

import numpy as np

from cosmology import PLANCK18, Cosmology
from framework import A0_SPARC_LOCAL, G_NEWTON
from btfr import M_SUN_KG, KM_PER_M, _E_eff


# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------
# kpc <-> m: same Mpc <-> km factor as in cosmology.py, scaled by 1e3 (kpc/Mpc).
M_PER_KPC = 3.085677581491367e19  # m per kpc  (= km per Mpc, by coincidence of unit factors)


# ---------------------------------------------------------------------------
# Anchored a_0(z)
# ---------------------------------------------------------------------------

def a0_of_z_anchored(z, cosmo: Cosmology = PLANCK18,
                     a0_local: float = A0_SPARC_LOCAL) -> float:
    """a_0(z) under the anchored convention a_0(0) = a_0_SPARC exactly.

    a_0(z) = a_0_SPARC * E(z)/E(0).
    """
    return a0_local * _E_eff(z, cosmo=cosmo)


# ---------------------------------------------------------------------------
# §5.1 MOND transition radius
# ---------------------------------------------------------------------------

def r_M(M_b_solar: float, z: float = 0.0,
        cosmo: Cosmology = PLANCK18,
        a0_local: float = A0_SPARC_LOCAL) -> float:
    """MOND transition radius in kpc, eq. (5.1):

        r_M(z) = sqrt( G * M_b / a_0(z) ).
    """
    M_b_kg = M_b_solar * M_SUN_KG
    a0_z = a0_of_z_anchored(z, cosmo=cosmo, a0_local=a0_local)
    r_m = np.sqrt(G_NEWTON * M_b_kg / a0_z)   # meters
    return r_m / M_PER_KPC                    # kpc


def r_M_ratio(z, cosmo: Cosmology = PLANCK18) -> float:
    """r_M(z) / r_M(0) = 1 / sqrt( E(z)/E(0) ), eq. (5.2) anchored."""
    return 1.0 / np.sqrt(_E_eff(z, cosmo=cosmo))


# ---------------------------------------------------------------------------
# §5 asymptotic flat velocity
# ---------------------------------------------------------------------------

def v_flat(M_b_solar: float, z: float = 0.0,
           cosmo: Cosmology = PLANCK18,
           a0_local: float = A0_SPARC_LOCAL) -> float:
    """Asymptotic flat velocity in km/s, deep-MOND BTFR (4.1):

        v_flat(z) = ( G * M_b * a_0(z) )^(1/4).
    """
    M_b_kg = M_b_solar * M_SUN_KG
    a0_z = a0_of_z_anchored(z, cosmo=cosmo, a0_local=a0_local)
    v_si = (G_NEWTON * M_b_kg * a0_z) ** 0.25
    return v_si * KM_PER_M


def v_flat_ratio(z, cosmo: Cosmology = PLANCK18) -> float:
    """v_flat(z) / v_flat(0) = ( E(z)/E(0) )^(1/4), eq. (5.3) anchored."""
    return _E_eff(z, cosmo=cosmo) ** 0.25


# ---------------------------------------------------------------------------
# §5.2 Table 3 archetypes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Archetype:
    name: str
    M_b: float       # solar masses
    R_d: float       # disk scale length, kpc (informational)


ARCHETYPES: List[Archetype] = [
    Archetype("Dwarf (DDO 154-like)",     3.5e8,  0.5),
    Archetype("Sub-L* (NGC 2403-like)",   1.0e10, 1.6),
    Archetype("L* (NGC 6946-like)",       6.0e10, 2.5),
    Archetype("Giant (UGC 2885-like)",    1.5e11, 8.0),
]


@dataclass(frozen=True)
class Table3Row:
    name: str
    M_b: float
    R_d: float
    rM0: float       # kpc
    rM2: float       # kpc
    v0: float        # km/s
    v2: float        # km/s


def table3_rows(cosmo: Cosmology = PLANCK18,
                a0_local: float = A0_SPARC_LOCAL) -> List[Table3Row]:
    rows: List[Table3Row] = []
    for arch in ARCHETYPES:
        rows.append(Table3Row(
            name=arch.name,
            M_b=arch.M_b,
            R_d=arch.R_d,
            rM0=r_M(arch.M_b, 0.0, cosmo=cosmo, a0_local=a0_local),
            rM2=r_M(arch.M_b, 2.0, cosmo=cosmo, a0_local=a0_local),
            v0=v_flat(arch.M_b, 0.0, cosmo=cosmo, a0_local=a0_local),
            v2=v_flat(arch.M_b, 2.0, cosmo=cosmo, a0_local=a0_local),
        ))
    return rows


# ---------------------------------------------------------------------------
# Verification driver
# ---------------------------------------------------------------------------

def main():
    all_match = True
    mismatches = []

    # ----- §5.1 fractional r_M shifts -----
    print("§5.1 fractional MOND-radius shifts (mass-independent):")
    fractional_targets = {
        2.0: 0.574,
        5.0: 0.347,
        10.0: 0.221,
    }
    for z, target in fractional_targets.items():
        val = float(r_M_ratio(z))
        rounded = round(val, 3)
        ok = rounded == target
        flag = "ok" if ok else "MISMATCH"
        print(f"  z = {z:5.1f}: r_M(z)/r_M(0) = {val:.6f}  "
              f"-> {rounded:.3f}  (paper: {target:.3f})  [{flag}]")
        if not ok:
            all_match = False
            mismatches.append((f"§5.1 r_M ratio z={z}", target, val, 3))
    print()

    # Caption shifts (r_M(2)/r_M(0) = 0.574 already covered; v_flat(2)/v_flat(0))
    v_ratio_2 = float(v_flat_ratio(2.0))
    v_ratio_2_rounded = round(v_ratio_2, 3)
    target_v_ratio = 1.320
    ok_v_ratio = v_ratio_2_rounded == target_v_ratio
    flag = "ok" if ok_v_ratio else "MISMATCH"
    print(f"§5.2 caption: v_flat(2)/v_flat(0) = {v_ratio_2:.6f} "
          f"-> {v_ratio_2_rounded:.3f}  (paper: {target_v_ratio:.3f})  [{flag}]")
    if not ok_v_ratio:
        all_match = False
        mismatches.append(("§5.2 caption v_flat(2)/v_flat(0)",
                           target_v_ratio, v_ratio_2, 3))
    print()

    # ----- §5.2 Table 3 archetypes -----
    print("§5.2 Table 3 archetypes (paper precisions: r_M to 2 dp, v to 0 dp):")
    paper_table = {
        "Dwarf (DDO 154-like)":   (0.64,  0.37,  49,  64),
        "Sub-L* (NGC 2403-like)": (3.41,  1.96, 112, 148),
        "L* (NGC 6946-like)":     (8.35,  4.79, 176, 232),
        "Giant (UGC 2885-like)":  (13.20, 7.58, 221, 292),
    }
    header = (f"{'Archetype':28s} {'M_b':>10s} {'R_d':>5s} "
              f"{'r_M(0)':>7s} {'r_M(2)':>7s} {'v(0)':>5s} {'v(2)':>5s}")
    print(header)
    rows = table3_rows()
    for row in rows:
        prM0, prM2, pv0, pv2 = paper_table[row.name]
        # Display precisions:
        #   r_M to 2 dp
        #   v_flat to nearest integer (paper shows no decimals)
        rrM0 = round(row.rM0, 2)
        rrM2 = round(row.rM2, 2)
        rv0 = int(round(row.v0))
        rv2 = int(round(row.v2))

        ok_rM0 = rrM0 == prM0
        ok_rM2 = rrM2 == prM2
        ok_v0 = rv0 == pv0
        ok_v2 = rv2 == pv2
        marks = "".join("." if ok else "X" for ok in
                        (ok_rM0, ok_rM2, ok_v0, ok_v2))
        print(f"{row.name:28s} {row.M_b:10.2e} {row.R_d:5.1f} "
              f"{row.rM0:7.4f} {row.rM2:7.4f} "
              f"{row.v0:5.1f} {row.v2:5.1f}   {marks}")

        for label, ok, comp, paper, dp in [
            (f"{row.name} r_M(0)", ok_rM0, row.rM0, prM0, 2),
            (f"{row.name} r_M(2)", ok_rM2, row.rM2, prM2, 2),
            (f"{row.name} v(0)",   ok_v0,  row.v0,  pv0,  0),
            (f"{row.name} v(2)",   ok_v2,  row.v2,  pv2,  0),
        ]:
            if not ok:
                all_match = False
                mismatches.append((label, paper, comp, dp))

    print()
    print("Paper Table 3 (for reference):")
    print(header)
    for arch in ARCHETYPES:
        prM0, prM2, pv0, pv2 = paper_table[arch.name]
        print(f"{arch.name:28s} {arch.M_b:10.2e} {arch.R_d:5.1f} "
              f"{prM0:7.2f} {prM2:7.2f} {pv0:5d} {pv2:5d}")
    print()

    # ----- §5.3 figure asymptotes (L*) -----
    Lstar = next(a for a in ARCHETYPES if a.name.startswith("L*"))
    v0_L = v_flat(Lstar.M_b, 0.0)
    v2_L = v_flat(Lstar.M_b, 2.0)
    print(f"§5.3 Figure 3 L* asymptotes:")
    print(f"  v_flat(0) = {v0_L:.4f} km/s  (paper: 176)  "
          f"match: {int(round(v0_L)) == 176}")
    print(f"  v_flat(2) = {v2_L:.4f} km/s  (paper: 232)  "
          f"match: {int(round(v2_L)) == 232}")

    print()
    print("=" * 60)
    if all_match:
        print("ALL §5 NUMERICAL CLAIMS MATCH PAPER EXACTLY (at displayed precision).")
    else:
        print("DISAGREEMENTS FOUND:")
        for label, paper, comp, dp in mismatches:
            print(f"  {label}: paper={paper}, computed={comp}  (dp={dp})")


if __name__ == "__main__":
    main()
