"""§4.1 Forward-model bias sweep for the Übler+ 2017 KMOS3D BTFR comparison.

This module implements the four-bias-model Monte Carlo forward analysis
sketched in §4.1 of the paper. The question it answers:

    Can a literature-plausible combination of measurement systematics applied
    to a framework-BTFR-obeying galaxy population reproduce the non-monotonic
    KMOS3D pattern Δb_obs = (-0.443 dex at z=0.9, -0.270 dex at z=2.3)?

The framework's prediction is strictly monotonic, Δb_MIT(z) = -log10(E(z)/E(0)),
giving Δb_MIT(0.9) = -0.227 and Δb_MIT(2.3) = -0.540. The observed pattern
crosses this prediction (more negative at z=0.9, less negative at z=2.3).
A trend-shape rescue would require systematics that push Δb in opposite
directions at the two redshifts.

Pipeline
--------
1. Generate mock galaxies obeying the framework BTFR at each z, using
   literature-realistic distributions of M_b (log-normal), R_d (size-mass
   relation), and σ_0 (Wisnioski+ 2015 medians).
2. Apply the Übler thick-disk correction
       v_circ²(r) = v_rot²(r) + 2 σ_0² × r/R_d
   at the characteristic measurement radius r = 2.2 R_d.
3. Apply the published selection cut v_rot,max / σ_0 > 4.4.
4. Sweep four bias models over their literature-plausible amplitude ranges:
     - RPU: radial-position uncertainty (recovered v_rot < true v_rot)
     - BS:  beam-smearing of σ_0 (apparent σ inflated)
     - NA:  non-asymptotic rotation curves (v_rot,max < v_flat)
     - SI:  selection-induced mass/threshold shift
   Each amplitude is in [0, 1] over the literature-plausible range; the
   z-dependence of each effect follows the standard scaling with angular
   resolution and pressure-support fraction.
5. For each amplitude tuple (α, β, γ, δ), recover the fixed-slope BTFR
   zero-point at z=0.9 and z=2.3, compute Δb against the local baseline.
6. Find the combination minimizing the joint residual to the observed
   (Δb_obs(0.9), Δb_obs(2.3)) — the "closest combined model".

Numerical conventions
---------------------
Random seed = 42 (per scripts/README.md). Mock population N = 5000 per z
gives ~0.005 dex statistical floor on the recovered b_const. Grid search uses
n=11 nodes per dimension (11^4 = 14,641 combinations) followed by a local
refinement step around the grid winner.

The local-baseline b(z=0) is computed from the closed-form deep-MOND relation
A_BTFR(0) = 1/(G a_0(0)) (≈ 62.78 M_sun/(km/s)^4, paper §3.1) so it is not
subject to the same Monte Carlo noise as the per-z recovered values.

This is the formal companion to ubler_sigma_tension.py: that module quantifies
the disagreement between Δb_obs and Δb_MIT under three error budgets; this
module asks whether bias models can dissolve the disagreement at all.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

import numpy as np

from cosmology import PLANCK18, Cosmology
from framework import A0_SPARC_LOCAL, G_NEWTON


# ---------------------------------------------------------------------------
# Constants and inputs
# ---------------------------------------------------------------------------

SEED = 42
N_GAL = 5000

M_SUN_KG = 1.98892e30        # kg
KPC_M = 3.085677581491367e19  # m
KMS_M = 1.0e3                 # m/s

# Wisnioski+ 2015 (KMOS3D) median intrinsic σ_0 (km/s).
# Used as the sample-mean σ at each z; per-galaxy σ_0 is drawn from a
# log-normal around this mean with 0.10 dex scatter (Wisnioski 2015 fig. 11).
SIGMA0_AT_Z_KMS: Dict[float, float] = {
    0.0: 25.0,
    0.9: 49.0,
    2.3: 70.0,
}

# Übler+ 2017 KMOS3D BTFR fixed-slope zero-point offsets relative to local Lelli.
# Δb_obs(0.9) = -0.443 (Übler+ 2017 Tab. 4; the paper §4.1 σ-tension table is
# consistent with -0.443 to two-decimal display -0.44; ubler_sigma_tension.py
# pins the third digit). Δb_obs(2.3) = -0.27.
DELTA_B_OBS: Dict[float, float] = {0.9: -0.443, 2.3: -0.27}

# Übler thick-disk correction: characteristic measurement radius in units of R_d.
R_MEAS_OVER_RD = 2.2  # standard KMOS3D / Übler convention

# Übler selection cut on raw rotation-to-dispersion ratio.
V_OVER_SIGMA_CUT = 4.4

# Mock M_b distribution: log10(M_b/M_sun) ~ Normal(MEAN_LOG_MB, SIGMA_LOG_MB),
# truncated to [9.5, 11.5] to match the Übler+ 2017 KMOS3D mass range
# (which spans roughly that decade after the v/σ > 4.4 selection).
MEAN_LOG_MB = 10.5
SIGMA_LOG_MB = 0.40


# ---------------------------------------------------------------------------
# Framework anchor and predictions
# ---------------------------------------------------------------------------

def E_eff(z: float, cosmo: Cosmology = PLANCK18) -> float:
    """E(z) anchored so E(0) = 1 exactly (matches btfr._E_eff)."""
    return float(cosmo.E(z)) / float(cosmo.E(0.0))


def a0_at(z: float) -> float:
    """Framework prediction a_0(z) = a_0(0) E(z) [m/s²]."""
    return A0_SPARC_LOCAL * E_eff(z)


def delta_b_MIT(z: float) -> float:
    """Δb_MIT(z) = -log10(E_eff(z))  [dex]."""
    return -math.log10(E_eff(z))


def b_local_baseline_solar_kms4() -> float:
    """Local b(z=0) = log10(A_BTFR(0)) with A in M_sun/(km/s)^4 (~ 1.798 dex).

    This is the closed-form deep-MOND zero-point used as the reference for
    Δb computation; it agrees with the local SPARC normalization in §3.1.
    """
    A_si = 1.0 / (G_NEWTON * A0_SPARC_LOCAL)        # kg / (m/s)^4
    A_solar = A_si * 1.0e12 / M_SUN_KG              # M_sun / (km/s)^4
    return math.log10(A_solar)


# ---------------------------------------------------------------------------
# Mock galaxy population
# ---------------------------------------------------------------------------

@dataclass
class Population:
    z: float
    log_Mb_solar: np.ndarray
    M_b_kg: np.ndarray
    v_circ_true: np.ndarray   # m/s
    v_rot_true: np.ndarray    # m/s
    sigma_0_true: np.ndarray  # m/s
    R_d_m: np.ndarray         # m
    R_meas_m: np.ndarray      # m


def generate_population(z: float, n: int = N_GAL,
                        rng: Optional[np.random.Generator] = None) -> Population:
    """Mock galaxies obeying framework BTFR at z under literature-realistic
    distributions of M_b, R_d, and σ_0."""
    if rng is None:
        rng = np.random.default_rng(SEED + int(round(z * 100)))

    # Baryonic mass: log-normal, truncated to [9.5, 11.5] dex
    log_Mb = rng.normal(MEAN_LOG_MB, SIGMA_LOG_MB, n)
    log_Mb = np.clip(log_Mb, 9.5, 11.5)
    M_b = (10.0 ** log_Mb) * M_SUN_KG  # kg

    # Framework v_circ from deep-MOND BTFR
    a0 = a0_at(z)
    v_circ_true = (G_NEWTON * M_b * a0) ** 0.25  # m/s

    # Disk scale length (van der Wel+ 2014 size-mass relation, broad form):
    #   log10(R_d / kpc) = 0.30 (log_Mb - 10) + 0.50 + N(0, 0.15)
    log_Rd_kpc = 0.30 * (log_Mb - 10.0) + 0.50 + rng.normal(0.0, 0.15, n)
    R_d = (10.0 ** log_Rd_kpc) * KPC_M  # m
    R_meas = R_MEAS_OVER_RD * R_d

    # Velocity dispersion (Wisnioski+ 2015 medians + 0.10 dex scatter)
    sigma_mean_ms = SIGMA0_AT_Z_KMS[z] * KMS_M
    log_sigma = np.log10(sigma_mean_ms) + rng.normal(0.0, 0.10, n)
    sigma_0 = 10.0 ** log_sigma  # m/s
    sigma_0 = np.clip(sigma_0, 5 * KMS_M, 200 * KMS_M)

    # Decompose v_circ_true into v_rot_true via the Übler thick-disk relation
    # at r = R_meas.
    pressure_term = 2.0 * sigma_0 ** 2 * (R_meas / R_d)  # = 2 σ² × R_MEAS_OVER_RD
    v_rot_sq = v_circ_true ** 2 - pressure_term

    # Galaxies for which the predicted v_circ cannot accommodate the drawn σ_0
    # are dropped (would imply pressure-dominated systems outside Übler's
    # rotation-supported sample). With v/σ > 4.4 selection downstream, this
    # only trims a small fraction at high z.
    valid = v_rot_sq > (5 * KMS_M) ** 2
    return Population(
        z=z,
        log_Mb_solar=log_Mb[valid],
        M_b_kg=M_b[valid],
        v_circ_true=v_circ_true[valid],
        v_rot_true=np.sqrt(v_rot_sq[valid]),
        sigma_0_true=sigma_0[valid],
        R_d_m=R_d[valid],
        R_meas_m=R_meas[valid],
    )


# ---------------------------------------------------------------------------
# Bias model
# ---------------------------------------------------------------------------

# Literature-plausible maximum bias amplitudes. Each is the multiplicative
# (or additive) shift at z = 0.9; the z = 2.3 strength scales by f_z(z, kind).
#
# Sources / rationale:
#   RPU: typical KMOS3D measurement-radius systematic ~5% on v_rot at z=0.9
#         (Genzel+ 2014, Übler+ 2017 Sec. 3.2); ~1.7× worse at z=2.3 because
#         the angular-diameter scale grows the apparent beam-to-disk ratio.
#   BS:  beam-smearing inflates apparent σ_0 by up to ~30% at the lowest
#         resolutions in KMOS3D (Wisnioski+ 2015 §3, Burkert+ 2016 fig. 4);
#         scales ~linearly with PSF/R_d_apparent.
#   NA:  v_rot,max can underestimate v_flat by up to ~7% when the rotation
#         curve is sampled before the asymptote (Genzel+ 2017 §2.2); the
#         offset grows with z because high-z disks are more often caught
#         in their rising portion.
#   SI:  the v/σ > 4.4 cut moves the survivor mean by a population-dependent
#         amount; we let the threshold itself shift by ±15% as a proxy for
#         heterogeneous selection across the two redshift bins.
BIAS_MAX = {
    'rpu': 0.05,   # max v_rot reduction from radial-position bias at z=0.9
    'bs':  0.30,   # max σ_0 inflation from beam-smearing at z=0.9
    'na':  0.07,   # max v_rot,max suppression from non-asymptotic curves at z=0.9
    'si':  0.15,   # max fractional shift in v/σ threshold (sign-controlled)
}


def f_z(z: float, kind: str) -> float:
    """Z-dependent strength multiplier for a given bias.

    z = 0.9 -> 1.0 by definition. z = 2.3 multipliers reflect
    literature-typical scalings of each systematic with redshift:
      - resolution-driven biases (rpu, bs, na) scale ~linearly to ~quadratically
        with the apparent angular size of the galaxy, ~ (1+z)^0.7 over this
        range -> a factor of ~1.6-1.9 from z=0.9 to z=2.3.
      - selection-driven biases (si) are not z-scaled here; we instead let
        the SI amplitude be signed to allow asymmetric impact across z.
    """
    if z == 0.9:
        return 1.0
    if z == 2.3:
        if kind == 'rpu':
            return 1.7
        if kind == 'bs':
            return 1.9
        if kind == 'na':
            return 1.7
        if kind == 'si':
            return 1.0
    raise ValueError(f"f_z table only defined at z in {{0.9, 2.3}}, got z={z}")


def apply_biases_and_recover(pop: Population,
                             alpha: float, beta: float,
                             gamma: float, delta: float,
                             si_sign: int = +1) -> float:
    """Apply (RPU, BS, NA, SI) biases at amplitudes (α, β, γ, δ) ∈ [0, 1] and
    recover the fixed-slope BTFR zero-point b_recov in M_sun/(km/s)^4 units.

    The si_sign argument lets the selection-cut shift go either direction
    (looser or stricter); the "closest combined" search optimizes both signs.

    Returns NaN if fewer than 50 galaxies survive the selection cut.
    """
    z = pop.z
    rpu = alpha * BIAS_MAX['rpu'] * f_z(z, 'rpu')
    bs = beta * BIAS_MAX['bs'] * f_z(z, 'bs')
    na = gamma * BIAS_MAX['na'] * f_z(z, 'na')
    si = si_sign * delta * BIAS_MAX['si'] * f_z(z, 'si')

    # RPU + NA both reduce recovered v_rot
    v_rot_obs = pop.v_rot_true * (1.0 - rpu) * (1.0 - na)

    # BS inflates apparent σ_0
    sigma_obs = pop.sigma_0_true * (1.0 + bs)

    # Übler thick-disk correction with observed σ at the measurement radius
    pressure_term_obs = 2.0 * sigma_obs ** 2 * (pop.R_meas_m / pop.R_d_m)
    v_circ_obs = np.sqrt(np.maximum(v_rot_obs ** 2 + pressure_term_obs,
                                    (5 * KMS_M) ** 2))

    # Selection cut on raw v_rot,obs / σ_obs, with SI shift on the threshold
    threshold = V_OVER_SIGMA_CUT * (1.0 + si)
    keep = v_rot_obs / sigma_obs > threshold

    if int(keep.sum()) < 50:
        return float('nan')

    log_Mb = pop.log_Mb_solar[keep]
    log_v_kms = np.log10(v_circ_obs[keep] / KMS_M)
    return float(np.median(log_Mb - 4.0 * log_v_kms))


# ---------------------------------------------------------------------------
# Sweep and closest-combined search
# ---------------------------------------------------------------------------

def recover_delta_b(pops: Dict[float, Population],
                    alpha: float, beta: float, gamma: float, delta: float,
                    si_sign: int) -> Optional[Dict[float, float]]:
    """Compute Δb_recov(z) at each z in pops for the given amplitudes."""
    b0 = b_local_baseline_solar_kms4()
    out: Dict[float, float] = {}
    for z, pop in pops.items():
        b = apply_biases_and_recover(pop, alpha, beta, gamma, delta,
                                     si_sign=si_sign)
        if math.isnan(b):
            return None
        out[z] = b - b0
    return out


def joint_residual(deltas: Dict[float, float]) -> float:
    """Sum of squared residuals between recovered and observed Δb."""
    return sum((deltas[z] - DELTA_B_OBS[z]) ** 2 for z in deltas)


def grid_search(pops: Dict[float, Population],
                n_grid: int = 11, verbose: bool = False,
                criterion: str = 'joint',
                z_focus: Optional[float] = None,
                direction: str = 'closer',
                ) -> Tuple[Tuple[float, float, float, float, int],
                           Dict[float, float], float]:
    """Coarse grid over (α, β, γ, δ) ∈ [0, 1]^4 × {-1, +1} for SI sign.

    criterion:
      'joint'        - minimize joint L2 residual to (Δb_obs(0.9), Δb_obs(2.3))
      'per_z_closer' - minimize |Δb_recov(z_focus) - Δb_obs(z_focus)|
                       (single-bin closest-fit; ignores the other bin)
      'per_z_extreme'- maximize bias-induced shift toward `direction` at z_focus
                       ('more_negative' or 'less_negative'); ignores other bin
    """
    grid = np.linspace(0.0, 1.0, n_grid)
    best_score = math.inf
    best_params: Tuple[float, float, float, float, int] = (0, 0, 0, 0, +1)
    best_deltas: Dict[float, float] = {}
    n_eval = 0
    n_skip = 0
    for si_sign in (+1, -1):
        for a in grid:
            for b in grid:
                for c in grid:
                    for d in grid:
                        deltas = recover_delta_b(pops, a, b, c, d, si_sign)
                        n_eval += 1
                        if deltas is None:
                            n_skip += 1
                            continue
                        if criterion == 'joint':
                            score = joint_residual(deltas)
                        elif criterion == 'per_z_closer':
                            assert z_focus is not None
                            score = abs(deltas[z_focus] - DELTA_B_OBS[z_focus])
                        elif criterion == 'per_z_extreme':
                            assert z_focus is not None
                            if direction == 'more_negative':
                                score = deltas[z_focus]  # smaller = more negative
                            else:
                                score = -deltas[z_focus]  # larger = less negative
                        else:
                            raise ValueError(f"unknown criterion {criterion!r}")
                        if score < best_score:
                            best_score = score
                            best_params = (float(a), float(b), float(c),
                                           float(d), int(si_sign))
                            best_deltas = deltas
    if verbose:
        print(f"  grid search ({criterion}): {n_eval} evals, "
              f"{n_skip} skipped (no survivors)")
    return best_params, best_deltas, best_score


def refine(pops: Dict[float, Population],
           best_params: Tuple[float, float, float, float, int],
           half_width: float = 0.10, n_grid: int = 11
           ) -> Tuple[Tuple[float, float, float, float, int],
                      Dict[float, float], float]:
    """Local refinement around the coarse-grid winner; SI sign held fixed."""
    a0, b0_, c0, d0, si_sign = best_params
    grid_a = np.linspace(max(0.0, a0 - half_width),
                         min(1.0, a0 + half_width), n_grid)
    grid_b = np.linspace(max(0.0, b0_ - half_width),
                         min(1.0, b0_ + half_width), n_grid)
    grid_c = np.linspace(max(0.0, c0 - half_width),
                         min(1.0, c0 + half_width), n_grid)
    grid_d = np.linspace(max(0.0, d0 - half_width),
                         min(1.0, d0 + half_width), n_grid)
    best_res = math.inf
    best_params_out = best_params
    best_deltas: Dict[float, float] = {}
    for a in grid_a:
        for b in grid_b:
            for c in grid_c:
                for d in grid_d:
                    deltas = recover_delta_b(pops, a, b, c, d, si_sign)
                    if deltas is None:
                        continue
                    r = joint_residual(deltas)
                    if r < best_res:
                        best_res = r
                        best_params_out = (float(a), float(b), float(c),
                                           float(d), int(si_sign))
                        best_deltas = deltas
    return best_params_out, best_deltas, best_res


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def baseline_check(pops: Dict[float, Population]) -> Dict[float, float]:
    """Sanity check: with zero bias, recovered Δb should equal Δb_MIT."""
    b0 = b_local_baseline_solar_kms4()
    out: Dict[float, float] = {}
    for z, pop in pops.items():
        b = apply_biases_and_recover(pop, 0.0, 0.0, 0.0, 0.0, si_sign=+1)
        out[z] = b - b0
    return out


# ---------------------------------------------------------------------------
# Canonical (committed) results for self-verification
# ---------------------------------------------------------------------------
# Generated by the deterministic seed=42, N_GAL=5000, n_grid=11 run. Used by
# verify.py to detect numerical drift; tolerances reflect the Monte Carlo
# floor (per-z medians of 5,000 mocks have ~0.005 dex noise).
EXPECTED = {
    'zero_bias': {0.9: -0.227, 2.3: -0.540},        # matches Δb_MIT exactly
    'joint_min': {0.9: -0.166, 2.3: -0.429},        # joint L2 minimum
    'joint_min_residual_dex': 0.319,
    # Single-bin closest-fit points: bias mix dragging Δb_recov(z) closest
    # to Δb_obs(z), regardless of the residual at the other bin.
    'best_z09_closer': {0.9: -0.292, 2.3: -0.648},  # 0.151 dex residual at 0.9
    'best_z23_closer': {0.9: -0.071, 2.3: -0.270},  # 0.000 dex residual at 2.3
    # Paper §4.1 narrative anchor: at the z=0.9-best params, the gap to
    # observed at z=0.9 is 0.151 dex; the simultaneous z=2.3 gap is -0.378.
    'paper_z09_residual_dex': +0.151,
    'paper_z23_residual_at_z09_best_dex': -0.378,
}
TOL_DEX = 0.020  # 20 mdex MC-floor + grid quantization tolerance


def _assert_close(label: str, computed: float, expected: float, tol: float):
    ok = abs(computed - expected) <= tol
    flag = "ok" if ok else "MISMATCH"
    print(f"  {label:<45s} computed={computed:+.3f}  "
          f"expected={expected:+.3f}  Δ={computed - expected:+.3f}  {flag}")
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n', 1)[0])
    parser.add_argument('--n-grid', type=int, default=11,
                        help='nodes per axis in the coarse grid (default 11)')
    parser.add_argument('--n-gal', type=int, default=N_GAL,
                        help='mock galaxies per redshift bin (default 5000)')
    parser.add_argument('--seed', type=int, default=SEED,
                        help='random seed (default 42)')
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args()

    z_list = (0.9, 2.3)
    print("§4.1 Übler+ 2017 KMOS3D BTFR forward-model bias sweep")
    print("=" * 72)
    print(f"Seed                 : {args.seed}")
    print(f"Mock galaxies / bin  : {args.n_gal}")
    print(f"Coarse grid (per dim): {args.n_grid}  ({args.n_grid}^4 × 2 signs"
          f" = {args.n_grid**4 * 2} evaluations)")
    print()

    # Generate populations once.
    rngs = {z: np.random.default_rng(args.seed + int(round(z * 100)))
            for z in z_list}
    pops = {z: generate_population(z, n=args.n_gal, rng=rngs[z])
            for z in z_list}
    for z, pop in pops.items():
        print(f"  z = {z}: {len(pop.log_Mb_solar)} mock galaxies after"
              f" σ-feasibility cut")
    print()

    # Predictions and observations
    print("Framework predictions and KMOS3D observations:")
    print(f"{'z':>5} {'Δb_MIT':>10} {'Δb_obs':>10} {'residual (obs-MIT)':>22}")
    for z in z_list:
        d_mit = delta_b_MIT(z)
        d_obs = DELTA_B_OBS[z]
        print(f"{z:5.2f} {d_mit:>+10.3f} {d_obs:>+10.3f}"
              f" {d_obs - d_mit:>+22.3f}")
    print()

    # Sanity: zero-bias should recover Δb_MIT(z) within Monte Carlo noise.
    sanity = baseline_check(pops)
    print("Zero-bias sanity check (should match Δb_MIT to within MC noise):")
    print(f"{'z':>5} {'Δb_recov(0,0,0,0)':>20} {'Δb_MIT':>10}"
          f" {'residual':>10}")
    for z in z_list:
        r = sanity[z] - delta_b_MIT(z)
        print(f"{z:5.2f} {sanity[z]:>+20.3f} {delta_b_MIT(z):>+10.3f}"
              f" {r:>+10.3f}")
    print()

    # ---- Joint L2 minimum --------------------------------------------------
    print("[1/3] Joint L2 minimum (closest combined model)")
    coarse_params, coarse_deltas, coarse_res = grid_search(
        pops, n_grid=args.n_grid, verbose=not args.quiet,
        criterion='joint')
    fine_params, fine_deltas, fine_res = refine(
        pops, coarse_params, half_width=1.0 / max(args.n_grid - 1, 1),
        n_grid=args.n_grid)
    if fine_res < coarse_res:
        joint_params, joint_deltas, joint_res = (
            fine_params, fine_deltas, fine_res)
    else:
        joint_params, joint_deltas, joint_res = (
            coarse_params, coarse_deltas, coarse_res)
    print(f"  params : (α={joint_params[0]:.3f}, β={joint_params[1]:.3f},"
          f" γ={joint_params[2]:.3f}, δ={joint_params[3]:.3f},"
          f" sign={joint_params[4]:+d})")
    for z in z_list:
        d_obs = DELTA_B_OBS[z]
        d_rec = joint_deltas[z]
        print(f"  z = {z}: Δb_recov = {d_rec:+.3f},"
              f" residual to obs = {d_rec - d_obs:+.3f} dex")
    print(f"  joint L2 distance to observed = "
          f"{math.sqrt(joint_res):.3f} dex")
    print()

    # ---- Per-bin "closest single bin" fits ---------------------------------
    print("[2/3] Best single-bin fit (cherry-pick z=0.9)")
    p09, d09, _ = grid_search(pops, n_grid=args.n_grid,
                              criterion='per_z_closer', z_focus=0.9,
                              verbose=not args.quiet)
    print(f"  params : (α={p09[0]:.3f}, β={p09[1]:.3f},"
          f" γ={p09[2]:.3f}, δ={p09[3]:.3f}, sign={p09[4]:+d})")
    for z in z_list:
        print(f"  z = {z}: Δb_recov = {d09[z]:+.3f},"
              f" residual to obs = {d09[z] - DELTA_B_OBS[z]:+.3f} dex")
    print()

    print("[3/3] Best single-bin fit (cherry-pick z=2.3)")
    p23, d23, _ = grid_search(pops, n_grid=args.n_grid,
                              criterion='per_z_closer', z_focus=2.3,
                              verbose=not args.quiet)
    print(f"  params : (α={p23[0]:.3f}, β={p23[1]:.3f},"
          f" γ={p23[2]:.3f}, δ={p23[3]:.3f}, sign={p23[4]:+d})")
    for z in z_list:
        print(f"  z = {z}: Δb_recov = {d23[z]:+.3f},"
              f" residual to obs = {d23[z] - DELTA_B_OBS[z]:+.3f} dex")
    print()

    # ---- Self-verification against committed canonical values -------------
    print("Self-verification against canonical seed=42 / N=5000 / n_grid=11"
          " values:")
    all_ok = True
    all_ok &= _assert_close("zero-bias  z=0.9 Δb",   sanity[0.9],
                            EXPECTED['zero_bias'][0.9], TOL_DEX)
    all_ok &= _assert_close("zero-bias  z=2.3 Δb",   sanity[2.3],
                            EXPECTED['zero_bias'][2.3], TOL_DEX)
    all_ok &= _assert_close("joint min  z=0.9 Δb",   joint_deltas[0.9],
                            EXPECTED['joint_min'][0.9], TOL_DEX)
    all_ok &= _assert_close("joint min  z=2.3 Δb",   joint_deltas[2.3],
                            EXPECTED['joint_min'][2.3], TOL_DEX)
    all_ok &= _assert_close("z=0.9-best Δb at z=0.9", d09[0.9],
                            EXPECTED['best_z09_closer'][0.9], TOL_DEX)
    all_ok &= _assert_close("z=2.3-best Δb at z=2.3", d23[2.3],
                            EXPECTED['best_z23_closer'][2.3], TOL_DEX)
    print()
    print("=" * 72)
    if all_ok:
        print("ALL §4.1 FORWARD-MODEL CLAIMS MATCH PAPER (within Monte Carlo "
              "tolerance).")
    else:
        print("DISAGREEMENT — see rows above.")
    print()
    print("Forward-model conclusion (paper §4.1):")
    r09 = joint_deltas[0.9] - DELTA_B_OBS[0.9]
    r23 = joint_deltas[2.3] - DELTA_B_OBS[2.3]
    print(f"  Across the literature-plausible bias-amplitude box, no")
    print(f"  combination reproduces the observed (-0.443, -0.270) pattern.")
    print(f"  At joint L2 minimum: Δb_recov = ({joint_deltas[0.9]:+.3f},"
          f" {joint_deltas[2.3]:+.3f}); residuals ({r09:+.3f},"
          f" {r23:+.3f}) dex.")
    print(f"  Even cherry-picking the z=0.9 bin alone, the bias models can")
    print(f"  shift Δb_recov(0.9) only to {d09[0.9]:+.3f} (gap to obs:"
          f" {d09[0.9] - DELTA_B_OBS[0.9]:+.3f} dex).")


if __name__ == '__main__':
    main()
