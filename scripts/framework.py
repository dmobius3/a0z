"""
framework.py
============

Bounded-topology framework primitives used throughout the paper:
the phase operator C(Theta), the edge-mode hierarchy normalization N_H(z),
the Fibonacci-well assignments at z = 0, and the scaling-law evaluation
that turns those into dimensional predictions.

References (paper section numbers are from a0-evolution-paper.md):
  - C(Theta) = 2 sin^2(pi Theta), with Theta on the 120-domain native to
    the binary-icosahedral quotient S^3 / 2I             (paper eq. (2.1))
  - Edge-mode scaling law:
        A_edge(z) / A_P = C(Theta_A) * N_H(z)            (eq. 2.1)
    with N_H(z) defined through the calibration relation
        H(z) t_P = C(34/120) * N_H(z)                    (eq. 2.3)
    so that
        N_H(z) = H(z) t_P / C(34/120).
  - Prediction for a_0(z):
        a_0(z) / a_P = C(13/120) * N_H(z)                (eq. 2.2)
    where a_P = c / t_P. Substituting (2.3) into (2.2) gives the
    structural ratio
        a_0(z) / (c H(z)) = C(13/120) / C(34/120) = 0.1845   (eq. 2.4)

Verification targets (per scripts/README.md):
  - C(13/120) / C(34/120) = 0.1845
  - Observed Milgrom ratio a_0 / (c H_0) = 0.1833
    (using a_0 = 1.20e-10 m/s^2 from SPARC and H_0 = 67.4 km/s/Mpc)

This module imports the cosmology only to evaluate H(z) for downstream
consumers; the §2 ratio itself is a pure topological number and does
not depend on cosmological parameters.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from cosmology import PLANCK18, Cosmology, KM_PER_MPC


# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
# CODATA 2018 values; only c, h-bar, G enter the Planck units used in §2.

C_LIGHT = 2.99792458e8                  # m / s, exact
HBAR = 1.054571817e-34                  # J s
G_NEWTON = 6.67430e-11                  # m^3 / (kg s^2)

# Planck units derived from (c, hbar, G):
#   t_P = sqrt(hbar G / c^5)
#   a_P = c / t_P  =  c * sqrt(c^5 / (hbar G))  =  sqrt(c^7 / (hbar G))
T_PLANCK = np.sqrt(HBAR * G_NEWTON / C_LIGHT**5)         # ~5.39e-44 s
A_PLANCK = C_LIGHT / T_PLANCK                            # ~5.56e51 m/s^2

# SPARC local acceleration scale, used only in the §2 absolute-value check.
A0_SPARC_LOCAL = 1.20e-10                                # m / s^2


# ---------------------------------------------------------------------------
# Phase operator on the 120-domain (§2)
# ---------------------------------------------------------------------------

DOMAIN_SIZE = 120  # binary-icosahedral quotient S^3 / 2I


def C_phase(theta):
    """
    Phase operator C(Theta) = 2 sin^2(pi Theta).

    `theta` may be a scalar or array of phase positions in [0, 1].
    Returns the unit-mean-normalized squared modulus of the anti-periodic
    Möbius ground mode (paper §2).
    """
    theta = np.asarray(theta, dtype=float)
    return 2.0 * np.sin(np.pi * theta) ** 2


def C_well(k: int) -> float:
    """C evaluated at the lattice point k / 120, for k in {0, ..., 119}."""
    if not (0 <= k < DOMAIN_SIZE):
        raise ValueError(f"well index k={k} outside 0..{DOMAIN_SIZE - 1}")
    return float(C_phase(k / DOMAIN_SIZE))


# ---------------------------------------------------------------------------
# Fibonacci-well assignments (Appendix A.3)
# ---------------------------------------------------------------------------
#
# The Fibonacci wells on the 120-domain are F_7..F_10 = {13, 21, 34, 55}/120.
# F_11 = 89 collapses to F_8 = 21 under the reflection symmetry
# C(k) = C(120 - k), so the four distinct Fibonacci wells give
# C(4, 2) = 6 unordered well-pairs (paper App A.3, §1, §2).
#
# Lambda sits at the antinode 60/120 (surface mode, n = 2). It is NOT a
# Fibonacci well; it is the topologically protected fixed point where
# d ln C / d Theta = 0 (paper App A.3, eligibility condition 3). Listed
# below in WELL_ASSIGNMENTS for completeness, but excluded from the
# Fibonacci subset.
#
# The two wells used by the present paper are:
#
#     a_0 -> 13/120     (edge mode, n = 1)
#     H   -> 34/120     (edge mode, n = 1)
#
# These are calibration inputs at z = 0, not derived; their structural
# ratio (2.4) is the testable claim of §2.

FIBONACCI_WELLS = (13, 21, 34, 55)

WELL_ASSIGNMENTS: Mapping[str, int] = {
    "a_0": 13,     # edge mode (Fibonacci well)
    "H":   34,     # edge mode (Fibonacci well)
    "Lambda": 60,  # surface mode antinode (App A.3) — not a Fibonacci well.
}

# Manifold-mode index n in the scaling law A/A_P = C(Theta) * N^n.
MODE_INDEX: Mapping[str, int] = {
    "a_0": 1,
    "H":   1,
    "Lambda": 2,
}


# ---------------------------------------------------------------------------
# Hierarchy normalization N_H(z) (eq. 2.3)
# ---------------------------------------------------------------------------

def _H_in_per_second(H_km_s_Mpc: float) -> float:
    """Convert H from km/s/Mpc to s^-1."""
    return H_km_s_Mpc / KM_PER_MPC      # (km/s) / km = 1/s


def N_H(z, cosmo: Cosmology = PLANCK18) -> float:
    """
    Edge-mode hierarchy normalization at epoch z, defined by the
    calibration relation (2.3):

        N_H(z) = H(z) * t_P / C(34/120).
    """
    H_si = _H_in_per_second(float(cosmo.H(z)))
    return H_si * T_PLANCK / C_well(34)


# ---------------------------------------------------------------------------
# Scaling-law evaluation (eq. 2.1, applied to a_0)
# ---------------------------------------------------------------------------

def a0_of_z(z, cosmo: Cosmology = PLANCK18) -> float:
    """
    Predicted MOND acceleration scale at epoch z, in m/s^2, from the
    edge-mode scaling law (2.2):

        a_0(z) = a_P * C(13/120) * N_H(z).

    Equivalently a_0(z) = (C(13/120)/C(34/120)) * c * H(z), which is the
    eq. (2.4) structural form.
    """
    return A_PLANCK * C_well(13) * N_H(z, cosmo=cosmo)


def a0_over_cH_predicted() -> float:
    """The eq. (2.4) structural ratio C(13/120) / C(34/120)."""
    return C_well(13) / C_well(34)


def a0_over_cH_observed(
    a0: float = A0_SPARC_LOCAL,
    H0_km_s_Mpc: float = PLANCK18.H0,
) -> float:
    """The observed Milgrom ratio a_0 / (c H_0), §2 absolute-value check."""
    H0_si = _H_in_per_second(H0_km_s_Mpc)
    return a0 / (C_LIGHT * H0_si)


# ---------------------------------------------------------------------------
# Self-check / verification
# ---------------------------------------------------------------------------

def main():
    # Phase-operator values at the two wells used in §2.
    C13 = C_well(13)
    C34 = C_well(34)
    ratio_pred = a0_over_cH_predicted()
    ratio_obs = a0_over_cH_observed()

    print("Phase-operator values at the §2 wells:")
    print(f"  C(13/120)         = {C13:.4f}   (paper: 0.2229)")
    print(f"  C(34/120)         = {C34:.4f}   (paper: 1.2079)")
    print()

    # Verification target 1: C(13/120) / C(34/120) = 0.1845
    target1 = 0.1845
    rounded1 = round(ratio_pred, 4)
    print(f"Verification target  C(13/120)/C(34/120) = {target1}")
    print(f"Computed             C(13/120)/C(34/120) = {ratio_pred:.10f}")
    print(f"Rounded to 4 dp                          = {rounded1}")
    print(f"Match (4 dp)         : {rounded1 == target1}")
    print()

    # Verification target 2: a_0 / (c H_0) = 0.1833
    target2 = 0.1833
    rounded2 = round(ratio_obs, 4)
    print(f"Verification target  a_0_obs / (c H_0_obs) = {target2}")
    print(f"  using a_0 = {A0_SPARC_LOCAL} m/s^2 (SPARC) and"
          f" H_0 = {PLANCK18.H0} km/s/Mpc (Planck 2018)")
    print(f"Computed             a_0 / (c H_0)         = {ratio_obs:.10f}")
    print(f"Rounded to 4 dp                            = {rounded2}")
    print(f"Match (4 dp)         : {rounded2 == target2}")
    print()

    # Cross-check: scaling-law evaluation reproduces the structural ratio
    # at z = 0 to round-off, independent of the chosen cosmology.
    a0_z0 = a0_of_z(0.0)
    H0_si = _H_in_per_second(PLANCK18.H0)
    cross = a0_z0 / (C_LIGHT * H0_si)
    print("Cross-check: scaling-law a_0(0) divided by c H_0:")
    print(f"  a_0(0) [predicted] = {a0_z0:.6e} m/s^2")
    print(f"  a_0(0) / (c H_0)   = {cross:.10f}  "
          f"(should equal {ratio_pred:.10f})")

    # Print the §2 percent-level disagreement between the structural
    # prediction and the observed ratio (paper quotes 0.8%).
    pct = 100.0 * (ratio_pred - ratio_obs) / ratio_obs
    print()
    print(f"Predicted vs observed Milgrom ratio: "
          f"{ratio_pred:.4f} vs {ratio_obs:.4f}  "
          f"({pct:+.2f}%)")


if __name__ == "__main__":
    main()
