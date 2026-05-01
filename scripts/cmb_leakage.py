"""
cmb_leakage.py
==============

Empirical bound on residual edge-mode leakage into cosmological perturbations
at recombination (paper §4.3).

Structural argument (§4.3)
--------------------------
The framework's selection rule (§2) assigns a_0 to the edge-mode sector
(n = 1, evolves through Omega_H(z)) and cosmological perturbations to the
space-mode sector (n = 3, governed by Omega_Lambda, no a_0(z) modification).
Under that assignment the framework's CMB prediction coincides with LambdaCDM's
and the structural leakage parameter is exactly epsilon = 0.

The §4.3 leakage bound asks how large any residual coupling could be without
violating Planck's first-peak amplitude precision.

Leakage ansatz (paper eq. in §4.3)
----------------------------------
    g_eff(R, z) = g_N(R, z) + epsilon * sqrt(g_N(R, z) * a_0(z)).

Under the Sachs-Wolfe relation Delta T / T ~ Phi / 3, a leakage at level
epsilon produces a fractional first-peak amplitude shift

    delta C_l^peak / C_l^peak  ~  epsilon * sqrt(a_0(z) / g_N).

Inverting:

    epsilon  <~  tolerance / sqrt(a_0(z) / g_N).

The most-constraining bound comes from the scale with the largest
sqrt(a_0/g_N), i.e. the smallest g_N — the sub-horizon BAO scale.

Paper inputs (§4.3)
-------------------
At z = 1090 (recombination):
  - a_0(z=1090) ~ 2.79e-6 m/s^2  (anchored a_0(z) = a_0_SPARC * E(z)/E(0))
  - sound-horizon perturbation:  R ~ 0.13 Mpc physical, g_N ~ 4.0e-11 m/s^2
  - sub-horizon BAO scale:       R ~ 0.05 Mpc physical, g_N ~ 1.5e-11 m/s^2

Planck first-peak amplitude precision (§4.3):
  - literal value      = 39 / 5733  ~  0.68%
  - conservative value = 0.5%   (the README headline tolerance)

Verification target (per scripts/README.md)
-------------------------------------------
  CMB leakage bound: epsilon  <~  1.2e-5 at 0.5% Planck tolerance.

In-paper §4.3 also quotes:
  - a_0(z=1090) ~ 2.79e-6 m/s^2
  - a_0(z=1090) / a_0(0)  ~ 23,000
  - sqrt(a_0/g_N) at sound-horizon scale = 264
  - sqrt(a_0/g_N) at sub-horizon BAO scale = 430
  - epsilon <~ 1.2e-5 at 0.5% tolerance
  - epsilon <~ 1.6e-5 at literal Planck precision (39/5733)

Imports
-------
  - cosmology.PLANCK18  -> E(z) at recombination
  - framework.A0_SPARC_LOCAL  -> 1.20e-10 m/s^2 (the anchored a_0 at z=0)

Notes on convention
-------------------
The paper uses the anchored convention a_0(z) = a_0_SPARC * E(z)/E(0)
throughout (§3.1, §3.1). The framework's structural form a_0(z) = a_P * C(13/120)
* N_H(z) (eq. 2.2) gives a_0(z=0) close to but slightly above the SPARC value;
the §§3-5 numerical predictions all use the SPARC-anchored form so that the
local calibration is the SPARC measurement. We follow the same convention
here, giving a_0(z=1090) = 2.79e-6 m/s^2 to three significant figures, the
value §4.3 reports.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from cosmology import PLANCK18, Cosmology
from framework import A0_SPARC_LOCAL


# ---------------------------------------------------------------------------
# Inputs (paper §4.3 narrative numbers)
# ---------------------------------------------------------------------------

Z_RECOMB = 1090.0  # redshift of recombination used by §4.3

# Newtonian gravitational acceleration of cosmological perturbations at the
# two scales §4.3 evaluates the bound on. SI units (m/s^2). These are the
# narrative values the paper text quotes; they enter the bound as ratios
# inside sqrt(a_0/g_N).
G_N_SOUND_HORIZON = 4.0e-11    # m/s^2   (R ~ 0.13 Mpc physical, sound horizon)
G_N_SUBHORIZON_BAO = 1.5e-11   # m/s^2   (R ~ 0.05 Mpc physical, sub-horizon BAO)

# Planck first-peak amplitude precision (§4.3).
# The "literal" value is the ratio §4.3 spells out: 39 / 5733.
# The "conservative" value is the rounded reading the README headline target
# uses (0.5%).
PLANCK_FIRST_PEAK_LITERAL_NUM = 39.0
PLANCK_FIRST_PEAK_LITERAL_DEN = 5733.0
PLANCK_FIRST_PEAK_LITERAL = (
    PLANCK_FIRST_PEAK_LITERAL_NUM / PLANCK_FIRST_PEAK_LITERAL_DEN
)
PLANCK_FIRST_PEAK_CONSERVATIVE = 0.005   # 0.5% conservative tolerance


# ---------------------------------------------------------------------------
# a_0 at recombination under the anchored convention
# ---------------------------------------------------------------------------

def a0_at_recombination(
    z: float = Z_RECOMB,
    cosmo: Cosmology = PLANCK18,
    a0_local: float = A0_SPARC_LOCAL,
) -> float:
    """
    Anchored a_0(z) = a_0_SPARC * E(z) / E(0)  (paper §3.1, applied here at
    z = 1090). Returns m/s^2.
    """
    Ez = float(cosmo.E(z))
    E0 = float(cosmo.E(0.0))
    return a0_local * Ez / E0


def a0_ratio_to_local(
    z: float = Z_RECOMB,
    cosmo: Cosmology = PLANCK18,
) -> float:
    """a_0(z) / a_0(0). Equals E(z) / E(0) under the anchored convention."""
    return float(cosmo.E(z)) / float(cosmo.E(0.0))


# ---------------------------------------------------------------------------
# Leakage geometry: sqrt(a_0(z) / g_N) at the two §4.3 scales
# ---------------------------------------------------------------------------

def sqrt_a0_over_gN(a0: float, gN: float) -> float:
    """sqrt(a_0(z) / g_N), the dimensionless prefactor in the §4.3 bound."""
    return float(np.sqrt(a0 / gN))


# ---------------------------------------------------------------------------
# §4.3 leakage bound
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LeakageBound:
    """Container for the §4.3 leakage-bound results at one tolerance."""
    tolerance: float            # fractional first-peak amplitude precision
    epsilon_sound_horizon: float
    epsilon_subhorizon_bao: float

    @property
    def epsilon_most_constraining(self) -> float:
        """The tighter (smaller) of the two scale-specific bounds."""
        return min(self.epsilon_sound_horizon, self.epsilon_subhorizon_bao)


def leakage_bound(
    tolerance: float,
    a0_rec: float,
    g_N_sound: float = G_N_SOUND_HORIZON,
    g_N_bao: float = G_N_SUBHORIZON_BAO,
) -> LeakageBound:
    """
    Invert the §4.3 ansatz at the given fractional tolerance.

        epsilon  <~  tolerance / sqrt(a_0(z=1090) / g_N).
    """
    eps_sound = tolerance / sqrt_a0_over_gN(a0_rec, g_N_sound)
    eps_bao = tolerance / sqrt_a0_over_gN(a0_rec, g_N_bao)
    return LeakageBound(
        tolerance=tolerance,
        epsilon_sound_horizon=eps_sound,
        epsilon_subhorizon_bao=eps_bao,
    )


# ---------------------------------------------------------------------------
# Self-check / verification against §4.3 quoted numbers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CheckRow:
    label: str
    paper: str
    computed: str
    matches: bool


def _format_two_sig(x: float) -> str:
    """Format x to two significant figures in scientific notation, like '1.2e-05'."""
    if x == 0:
        return "0"
    exp = int(np.floor(np.log10(abs(x))))
    mantissa = x / (10.0**exp)
    return f"{round(mantissa, 1):.1f}e{exp:+03d}"


def run_verification() -> Tuple[list, bool]:
    """Reproduce every §4.3 numerical claim and check it against the paper."""

    # 1. a_0 at recombination under the anchored convention.
    a0_rec = a0_at_recombination()
    a0_rec_paper = 2.79e-6
    rounded_a0 = round(a0_rec / 1e-6, 2) * 1e-6  # round to 3 sig figs in 1e-6
    a0_match = abs(rounded_a0 - a0_rec_paper) < 1e-9

    # 2. Ratio a_0(z=1090) / a_0(0)
    ratio = a0_ratio_to_local()
    # Paper says "approximately 23,000 times the local value" - 2 sig figs.
    ratio_rounded_2sf = round(ratio, -3)  # nearest 1000
    ratio_match = ratio_rounded_2sf == 23000.0

    # 3. sqrt(a_0/g_N) at sound-horizon scale
    s_sound = sqrt_a0_over_gN(a0_rec, G_N_SOUND_HORIZON)
    s_sound_match = round(s_sound) == 264

    # 4. sqrt(a_0/g_N) at sub-horizon BAO scale
    s_bao = sqrt_a0_over_gN(a0_rec, G_N_SUBHORIZON_BAO)
    # Paper rounds 431 -> 430 at 2 significant figures
    s_bao_match = round(s_bao, -1) == 430.0

    # 5. Literal Planck precision 39/5733 ~ 0.68%
    literal = PLANCK_FIRST_PEAK_LITERAL
    literal_match = round(literal * 100, 2) == 0.68

    # 6. Bound at 0.5% conservative tolerance
    bound_conservative = leakage_bound(PLANCK_FIRST_PEAK_CONSERVATIVE, a0_rec)
    eps_cons = bound_conservative.epsilon_most_constraining
    eps_cons_str = _format_two_sig(eps_cons)
    eps_cons_match = eps_cons_str == "1.2e-05"

    # 7. Bound at literal Planck precision 39/5733 ~ 0.68%
    bound_literal = leakage_bound(literal, a0_rec)
    eps_lit = bound_literal.epsilon_most_constraining
    eps_lit_str = _format_two_sig(eps_lit)
    eps_lit_match = eps_lit_str == "1.6e-05"

    # 8. Most-constraining scale: sub-horizon BAO (largest sqrt(a_0/g_N))
    most_constraining_is_bao = (
        bound_conservative.epsilon_subhorizon_bao
        < bound_conservative.epsilon_sound_horizon
    )

    rows = [
        CheckRow(
            "a_0(z=1090) [m/s^2]",
            f"~ 2.79e-6",
            f"{a0_rec:.4e}",
            a0_match,
        ),
        CheckRow(
            "a_0(z=1090) / a_0(0)",
            "~ 23,000",
            f"{ratio:.1f}",
            ratio_match,
        ),
        CheckRow(
            "sqrt(a_0/g_N) sound-horizon (g_N=4e-11, R~0.13 Mpc)",
            "264",
            f"{s_sound:.2f}",
            s_sound_match,
        ),
        CheckRow(
            "sqrt(a_0/g_N) sub-horizon BAO  (g_N=1.5e-11, R~0.05 Mpc)",
            "430",
            f"{s_bao:.2f}",
            s_bao_match,
        ),
        CheckRow(
            "Planck first-peak precision 39/5733",
            "~ 0.68%",
            f"{literal*100:.4f}%",
            literal_match,
        ),
        CheckRow(
            "epsilon  <=  1.2e-5  at 0.5% tolerance",
            "1.2e-05",
            eps_cons_str,
            eps_cons_match,
        ),
        CheckRow(
            "epsilon  <=  1.6e-5  at literal Planck precision",
            "1.6e-05",
            eps_lit_str,
            eps_lit_match,
        ),
        CheckRow(
            "Most-constraining scale is sub-horizon BAO",
            "True",
            f"{most_constraining_is_bao}",
            most_constraining_is_bao,
        ),
    ]
    all_match = all(r.matches for r in rows)
    return rows, all_match


def main():
    print("Module 8: cmb_leakage.py  (paper §4.3 leakage bound)")
    print("=" * 72)
    print()

    a0_rec = a0_at_recombination()
    print(f"Recombination redshift:           z = {Z_RECOMB:.0f}")
    print(f"E({Z_RECOMB:.0f})                          = {float(PLANCK18.E(Z_RECOMB)):.4f}")
    print(f"a_0(z=0) [SPARC anchor]           = {A0_SPARC_LOCAL:.3e} m/s^2")
    print(f"a_0(z=1090) anchored              = {a0_rec:.4e} m/s^2")
    print(f"a_0(z=1090) / a_0(z=0)            = {a0_ratio_to_local():.1f}")
    print()

    print("Per-scale leakage geometry (§4.3):")
    print(f"  sound-horizon   g_N = {G_N_SOUND_HORIZON:.2e} m/s^2"
          f"   sqrt(a_0/g_N) = {sqrt_a0_over_gN(a0_rec, G_N_SOUND_HORIZON):.2f}")
    print(f"  sub-horizon BAO g_N = {G_N_SUBHORIZON_BAO:.2e} m/s^2"
          f"   sqrt(a_0/g_N) = {sqrt_a0_over_gN(a0_rec, G_N_SUBHORIZON_BAO):.2f}")
    print()

    print("Planck first-peak amplitude precision:")
    print(f"  literal     = 39 / 5733 = {PLANCK_FIRST_PEAK_LITERAL*100:.4f}%"
          f"   (rounds to ~ 0.68%)")
    print(f"  conservative                  = {PLANCK_FIRST_PEAK_CONSERVATIVE*100:.4f}%"
          f"   (README headline tolerance)")
    print()

    bound_cons = leakage_bound(PLANCK_FIRST_PEAK_CONSERVATIVE, a0_rec)
    bound_lit = leakage_bound(PLANCK_FIRST_PEAK_LITERAL, a0_rec)

    print("Leakage bounds  (epsilon  <~  tolerance / sqrt(a_0/g_N)):")
    print(f"  at 0.5% tolerance:")
    print(f"     sound-horizon   eps  <=  {bound_cons.epsilon_sound_horizon:.4e}")
    print(f"     sub-horizon BAO eps  <=  {bound_cons.epsilon_subhorizon_bao:.4e}"
          f"   <- most constraining")
    print(f"  at literal Planck precision (39/5733):")
    print(f"     sound-horizon   eps  <=  {bound_lit.epsilon_sound_horizon:.4e}")
    print(f"     sub-horizon BAO eps  <=  {bound_lit.epsilon_subhorizon_bao:.4e}"
          f"   <- most constraining")
    print()

    print("Verification against §4.3 numerical claims:")
    print("-" * 72)
    rows, all_match = run_verification()
    for r in rows:
        flag = "PASS" if r.matches else "FAIL"
        print(f"  [{flag}]  {r.label}")
        print(f"           paper    : {r.paper}")
        print(f"           computed : {r.computed}")
    print("-" * 72)
    print(f"All §4.3 numerical claims match: {all_match}")
    print()

    print("Framework structural prediction:  epsilon = 0  exactly,")
    print("under the §2 selection rule's n = 3 space-mode assignment for")
    print("cosmological perturbations.  The empirical bound above shows the")
    print("structural prediction satisfies Planck consistency by every margin")
    print("Planck currently provides.")


if __name__ == "__main__":
    main()
