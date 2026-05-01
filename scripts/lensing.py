"""
lensing.py
==========

Galaxy-galaxy weak lensing predictions across cosmic time (paper §3.3).

Framework prediction (§3.3, eq. (3.3)):

    M_dyn(R, z) / M_b = (R / r_M(0)) * sqrt( E(z)/E(0) ),

derived by combining the Newtonian inversion (§3.3 proxy)

    M_dyn(R) = v_flat^2 R / G

with the deep-MOND BTFR (§3.1, eq. (3.1))

    v_flat^4 = G * M_b * a_0(z),

and the framework's anchored a_0(z) = a_0(0) * E(z)/E(0) from §3.

Conventions inherited from Modules 3 and 4:
  - Anchored a_0(z), so that a_0(0) = a_0_SPARC = 1.20e-10 m/s^2 exactly.
    Equivalently, r_M(0) and v_flat(0) reduce to the SPARC-calibrated
    archetype values from §3.2; the lensing fractional shift is
    sqrt( E(z)/E(0) ) under btfr._E_eff.

LambdaCDM comparison (Appendix B, also §3.3):
  Reference (M_halo, c) values for the L* archetype:

      z        M_halo [M_sun]    c
      0        1.5e12            7.5
      1        1.0e12            5.0
      2        7.0e11            3.5

  Closed-form NFW enclosed mass (Appendix B):

      M_NFW(<R) = M_halo * [ ln(1 + R/r_s) - (R/r_s)/(1 + R/r_s) ]
                          / [ ln(1 + c)     - c/(1 + c)         ],

  with scale radius r_s = R_vir/c and virial radius defined under the
  R_200 (200 * rho_crit(z)) convention:

      R_vir^3 = 3 * M_halo / (4 pi * 200 * rho_crit(0) * E(z)^2).

  The dynamical-mass-to-baryonic-mass ratio under Newtonian inversion is

      M_dyn(R, z) / M_b = (M_NFW(<R) + M_b) / M_b.

  Note: this script computes rho_crit(0) = 3 H_0^2 / (8 pi G) from the
  Module 1 H_0 (in km/s/Mpc) for full self-consistency. Appendix B quotes
  rho_crit(0) = 8.55e-27 kg/m^3 alongside H_0 = 67.4 km/s/Mpc; the
  computed value is 8.533e-27 kg/m^3, a 0.2% difference that is well
  below the §3.3 discriminator scale and shifts the Lambda-CDM ratios at
  the third decimal only.

  The un-normalized cosmology.E(z) is used for the LambdaCDM machinery,
  per the convention reminder in scripts/README.md: H(z) and rho_crit(z)
  use the raw Friedmann E(z), while the framework's a_0(z) and the
  sqrt(E(z)) lensing factor use the anchored E(z)/E(0).

§3.3 + Appendix B paper claims this script verifies:
  - §3.3 sqrt(E(z)) table:  1.150 / 1.338 / 1.741 at z = 0.5 / 1 / 2
    and percent excesses 15.0% / 33.8% / 74.1%.
  - §3.3 prose: LambdaCDM L* M_dyn/M_b at R = 100 kpc =
    14.0 / 13.6 / 13.8 at z = 0 / 1 / 2.
  - §3.3 prose: framework L* M_dyn/M_b at R = 100 kpc =
    11.98 / 16.03 / 20.86 at z = 0 / 1 / 2.
  - §3.3 prose: framework-vs-LambdaCDM factor 1.52 at z = 2.
  - §3.3 Table 4: framework L* predictions at R = 30 / 100 / 300 kpc and
    z = 0 / 0.5 / 1 / 2 / 5; values within each row scale as sqrt(E(z)).
  - §3.3 prose: 92% inferred dark-mass fraction at z = 0; 95% at z = 2,
    using R = 100 kpc framework values 11.98 and 20.86.
  - §3.3 Figure 5: vertical offset between curves equals
    log10(sqrt(E(z))) at each redshift.

Imports:
  - cosmology: PLANCK18, raw E(z), H_0, KM_PER_MPC
  - framework: G_NEWTON, A0_SPARC_LOCAL
  - btfr:      M_SUN_KG, _E_eff (anchored E(z)/E(0))
  - rotation_curves: r_M (kpc), M_PER_KPC

Coding standard (per scripts/README.md):
  Closed-form NFW M(<R); no halo-mass-function fitting and no quadrature.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Mapping, Tuple

import numpy as np

from cosmology import PLANCK18, Cosmology, KM_PER_MPC
from framework import A0_SPARC_LOCAL, G_NEWTON
from btfr import M_SUN_KG, _E_eff
from rotation_curves import r_M, M_PER_KPC


# ---------------------------------------------------------------------------
# Cosmology auxiliaries: critical density at z = 0 and z
# ---------------------------------------------------------------------------

def rho_crit_0_si(cosmo: Cosmology = PLANCK18) -> float:
    """Critical density at z = 0 in SI units (kg/m^3) from Module 1 H_0:
        rho_crit(0) = 3 H_0^2 / (8 pi G).
    """
    H0_si = cosmo.H0 / KM_PER_MPC          # s^-1
    return 3.0 * H0_si ** 2 / (8.0 * np.pi * G_NEWTON)


def rho_crit_z_si(z, cosmo: Cosmology = PLANCK18) -> float:
    """rho_crit(z) = rho_crit(0) * E(z)^2 with the un-normalized
    cosmology.E(z) (Friedmann), per the Appendix B prescription."""
    return rho_crit_0_si(cosmo=cosmo) * float(cosmo.E(z)) ** 2


# ---------------------------------------------------------------------------
# §3.3 Newtonian-inversion proxy and §3.3 framework prediction (eq. (3.3))
# ---------------------------------------------------------------------------

def M_dyn_over_Mb_framework(R_kpc: float, z: float, M_b_solar: float,
                            cosmo: Cosmology = PLANCK18,
                            a0_local: float = A0_SPARC_LOCAL) -> float:
    """Framework prediction (eq. (3.3)):

        M_dyn(R, z) / M_b  =  R / r_M(0) * sqrt( E(z)/E(0) ).

    Independent of M_b in the deep-MOND limit; the M_b dependence enters
    only through r_M(0) = sqrt(G M_b / a_0(0))."""
    rM0_kpc = r_M(M_b_solar, 0.0, cosmo=cosmo, a0_local=a0_local)
    return (R_kpc / rM0_kpc) * np.sqrt(_E_eff(z, cosmo=cosmo))


def fractional_lensing_enhancement(z, cosmo: Cosmology = PLANCK18) -> float:
    """sqrt( E(z)/E(0) ): the §3.3 universal lensing enhancement factor."""
    return float(np.sqrt(_E_eff(z, cosmo=cosmo)))


# ---------------------------------------------------------------------------
# Appendix B: closed-form NFW M(<R) + R_200 virial radius
# ---------------------------------------------------------------------------

def virial_radius_200c(M_halo_solar: float, z: float,
                       cosmo: Cosmology = PLANCK18) -> float:
    """R_200c in meters, from Appendix B definition

        R_vir^3 = 3 M_halo / (4 pi * 200 * rho_crit(0) * E(z)^2).

    Uses raw cosmology.E(z) (Module 1) and rho_crit(0) computed from H_0."""
    M_halo_kg = M_halo_solar * M_SUN_KG
    rho_z = rho_crit_z_si(z, cosmo=cosmo)        # = rho_crit(0)*E(z)^2
    R_vir_3 = 3.0 * M_halo_kg / (4.0 * np.pi * 200.0 * rho_z)
    return R_vir_3 ** (1.0 / 3.0)                # meters


def M_NFW_enclosed(R_m: float, M_halo_solar: float, c: float,
                   z: float, cosmo: Cosmology = PLANCK18) -> float:
    """Closed-form NFW M(<R) in kg, Appendix B eq.:

        M_NFW(<R) = M_halo * [ ln(1+x) - x/(1+x) ]
                            / [ ln(1+c) - c/(1+c) ],

    with x = R/r_s, r_s = R_vir/c. R_m is the aperture in meters."""
    M_halo_kg = M_halo_solar * M_SUN_KG
    R_vir = virial_radius_200c(M_halo_solar, z, cosmo=cosmo)
    r_s = R_vir / c
    x = R_m / r_s
    num = np.log(1.0 + x) - x / (1.0 + x)
    den = np.log(1.0 + c) - c / (1.0 + c)
    return M_halo_kg * num / den


def M_dyn_over_Mb_LambdaCDM(R_kpc: float, z: float,
                             M_b_solar: float,
                             M_halo_solar: float, c: float,
                             cosmo: Cosmology = PLANCK18) -> float:
    """LambdaCDM Newtonian-inversion proxy under Appendix B prescription:

        M_dyn(R, z) / M_b  =  ( M_NFW(<R; M_halo, c, z) + M_b ) / M_b.

    R_kpc is the physical aperture in kpc."""
    R_m = R_kpc * M_PER_KPC
    M_NFW = M_NFW_enclosed(R_m, M_halo_solar, c, z, cosmo=cosmo)
    M_b_kg = M_b_solar * M_SUN_KG
    return (M_NFW + M_b_kg) / M_b_kg


# ---------------------------------------------------------------------------
# Archetype baryonic masses and Appendix B (M_halo, c) for the L* archetype
# ---------------------------------------------------------------------------

# §3.3 / Table 4 / Appendix B: L* anchor.
M_B_LSTAR = 6.0e10  # M_sun

# Appendix B table: Moster-style SHMR + Duffy-style concentration at L*.
LSTAR_HALO_TABLE: Mapping[float, Tuple[float, float]] = {
    0.0: (1.5e12, 7.5),
    1.0: (1.0e12, 5.0),
    2.0: (7.0e11, 3.5),
}


# ---------------------------------------------------------------------------
# Verification driver
# ---------------------------------------------------------------------------

def main():
    all_match = True
    mismatches: List[Tuple[str, object, object]] = []

    def check(label: str, paper_value, computed_value, dp: int):
        nonlocal all_match
        rounded = round(computed_value, dp)
        if isinstance(paper_value, float):
            ok = rounded == paper_value
        else:
            ok = rounded == paper_value
        flag = "ok" if ok else "MISMATCH"
        print(f"  {label}: paper={paper_value}, computed={computed_value:.6f}, "
              f"rounded={rounded}  [{flag}]")
        if not ok:
            all_match = False
            mismatches.append((label, paper_value, computed_value))

    print("§3.3 universal lensing enhancement sqrt(E(z)):")
    sqE_targets = {0.5: 1.150, 1.0: 1.338, 2.0: 1.741}
    for z, target in sqE_targets.items():
        val = fractional_lensing_enhancement(z)
        check(f"sqrt(E(z)) at z={z}", target, val, 3)

    print()
    print("§3.3 percent excesses in M_dyn/M_b (== 100*(sqrt(E(z))-1)):")
    pct_targets = {0.5: 15.0, 1.0: 33.8, 2.0: 74.1}
    for z, target in pct_targets.items():
        val = (fractional_lensing_enhancement(z) - 1.0) * 100.0
        check(f"% excess at z={z}", target, val, 1)

    print()
    print("§3.3 / §3.3 framework L* M_dyn/M_b at R = 100 kpc:")
    framework_R100_targets = {0.0: 11.98, 1.0: 16.03, 2.0: 20.86}
    for z, target in framework_R100_targets.items():
        val = M_dyn_over_Mb_framework(100.0, z, M_B_LSTAR)
        check(f"framework L* @ R=100 kpc, z={z}", target, val, 2)

    print()
    print("§3.3 LambdaCDM L* M_dyn/M_b at R = 100 kpc (Appendix B):")
    lcdm_R100_targets = {0.0: 14.0, 1.0: 13.6, 2.0: 13.8}
    for z, target in lcdm_R100_targets.items():
        M_h, c = LSTAR_HALO_TABLE[z]
        val = M_dyn_over_Mb_LambdaCDM(100.0, z, M_B_LSTAR, M_h, c)
        check(f"LambdaCDM L* @ R=100 kpc, z={z}", target, val, 1)

    print()
    print("§3.3 framework-vs-LambdaCDM factor at z = 2 (R = 100 kpc):")
    M_h, c = LSTAR_HALO_TABLE[2.0]
    fw2 = M_dyn_over_Mb_framework(100.0, 2.0, M_B_LSTAR)
    lc2 = M_dyn_over_Mb_LambdaCDM(100.0, 2.0, M_B_LSTAR, M_h, c)
    check("framework / LambdaCDM at z=2", 1.52, fw2 / lc2, 2)

    print()
    print("§3.3 Table 4 (framework L* archetype, R = 30 / 100 / 300 kpc):")
    paper_table4 = {
        # (R_kpc, z): paper value at 2 dp
        (30,  0.0): 3.59,  (30,  0.5): 4.13,  (30,  1.0): 4.81,
        (30,  2.0): 6.26,  (30,  5.0): 10.35,
        (100, 0.0): 11.98, (100, 0.5): 13.77, (100, 1.0): 16.03,
        (100, 2.0): 20.86, (100, 5.0): 34.50,
        (300, 0.0): 35.93, (300, 0.5): 41.32, (300, 1.0): 48.08,
        (300, 2.0): 62.57, (300, 5.0): 103.50,
    }
    R_list = [30, 100, 300]
    z_list = [0.0, 0.5, 1.0, 2.0, 5.0]
    header_label = "R/z"
    print(f"  {header_label:>6s}", *[f"{z:>8.1f}" for z in z_list])
    for R in R_list:
        cells = []
        for z in z_list:
            val = M_dyn_over_Mb_framework(float(R), z, M_B_LSTAR)
            cells.append(val)
        print(f"  {R:>6d}", *[f"{v:>8.4f}" for v in cells])
        for z, val in zip(z_list, cells):
            target = paper_table4[(R, z)]
            check(f"Table4 R={R} kpc, z={z}", target, val, 2)

    print()
    print("§3.3 inferred dark-mass fraction at z = 0 and z = 2 "
          "(R = 100 kpc, L*):")
    fw0 = M_dyn_over_Mb_framework(100.0, 0.0, M_B_LSTAR)
    fw2 = M_dyn_over_Mb_framework(100.0, 2.0, M_B_LSTAR)
    dark_pct_0 = 100.0 * (1.0 - 1.0 / fw0)
    dark_pct_2 = 100.0 * (1.0 - 1.0 / fw2)
    check("dark-mass fraction at z=0 (paper: 92%)", 92, dark_pct_0, 0)
    check("dark-mass fraction at z=2 (paper: 95%)", 95, dark_pct_2, 0)

    print()
    print("§3.3 Figure 5 vertical offset = log10(sqrt(E(z))):")
    for z, sqE_target in [(0.0, 1.0), (1.0, 1.338), (2.0, 1.741)]:
        offset = np.log10(fractional_lensing_enhancement(z))
        target_offset = round(np.log10(sqE_target), 4)
        # Paper does not quote the offset numerically; just sanity check
        # that the framework curves are equally spaced in log10.
        print(f"  z={z}: log10(sqrt(E(z))) = {offset:.6f} "
              f"(target ~= {target_offset})")

    # Per scripts/README.md headline: ΛCDM L* M_dyn/M_b at R = 100 kpc:
    # 14.0, 13.6, 13.8 at z = 0, 1, 2.
    print()
    print("=" * 70)
    if all_match:
        print("ALL §3.3 + APPENDIX B NUMERICAL CLAIMS MATCH PAPER EXACTLY "
              "(at displayed precision).")
    else:
        print(f"DISAGREEMENTS FOUND ({len(mismatches)}):")
        for label, paper, comp in mismatches:
            print(f"  {label}: paper={paper}, computed={comp:.6f}")


if __name__ == "__main__":
    main()
