"""
make_figures.py
===============

Module 11 of the a0z analysis pipeline: generate Figures 1-5 of the paper
from the public APIs of the prior modules. No physics is rederived here;
this script is purely a presentation layer that calls into:

    cosmology       (Module 1) for E(z)
    framework       (Module 2) for a_0(z), C(Theta), etc.
    btfr            (Module 3) for A_BTFR(z) and the anchored E(z)/E(0)
    rotation_curves (Module 4) for r_M(z), v_flat(z), archetypes
    lensing         (Module 5) for M_dyn/M_b framework + ΛCDM, LSTAR_HALO_TABLE
    jwst_speedup    (Module 6) for the anchored E_eff helper (not used here)

Figures generated (paper sections in parentheses):

    Figure 1 (§3.2)  a_0(z)/a_0(0) = E(z)/E(0) from z=0 to z=15, log-y.
    Figure 2 (§4.3)  BTFR loci log10(M_b) vs log10(v_flat) at z = 0, 1, 2,
                     v_flat = 30..400 km/s, slope-4 lines anchored to
                     A_BTFR(0) = 62.78 M_sun/(km/s)^4.
    Figure 3 (§5.3)  L* archetype rotation curve at z = 0 and z = 2,
                     simple-MOND mu(x)=x/(1+x), point-mass g_N = G M_b / r^2,
                     plotted for r >= R_d (= 2.5 kpc).
    Figure 4 (§5.4)  r_M(z) for the four §5 archetypes, z = 0..5.
    Figure 5 (§6.4)  L* M_dyn(R)/M_b for R = 10..500 kpc at z = 0, 1, 2;
                     framework solid, ΛCDM (Appendix B halo table) dashed.

Each figure is saved to ../figures/figureN.{pdf,png}.

Each plot carries an in-script assertion against a paper-quoted value:
    Figure 1:  E(2)/E(0) = 3.033  (paper Table 1)
    Figure 2:  z=2 line passes through (log10(176)+0.32/4, log10(6e10))
               equivalently, log10(A_BTFR(2)) = log10(62.78/E(2)) = 1.316.
    Figure 3:  L* z=0 plateau -> 176 km/s; z=2 plateau -> 232 km/s.
    Figure 4:  L* r_M(0) = 8.35 kpc, r_M(2) = 4.79 kpc.
    Figure 5:  L* z=0 R=100 kpc framework value = 11.98 (paper §6.3).

Run from anywhere; output paths are resolved relative to the script's
parent directory.

Usage:
    python3 scripts/make_figures.py
"""

from __future__ import annotations

import os
import sys
from typing import Iterable, Tuple

import numpy as np
import matplotlib

# Use a non-interactive backend so the script runs in headless contexts.
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Make the prior modules importable when running as a script.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from cosmology import PLANCK18                                  # Module 1
from framework import A0_SPARC_LOCAL, G_NEWTON                  # Module 2
from btfr import (                                              # Module 3
    A_BTFR_solar_per_kms4,
    A_ratio,
    M_SUN_KG,
    KM_PER_M,
    _E_eff,
)
from rotation_curves import (                                   # Module 4
    ARCHETYPES,
    M_PER_KPC,
    a0_of_z_anchored,
    r_M,
    v_flat,
)
from lensing import (                                           # Module 5
    LSTAR_HALO_TABLE,
    M_B_LSTAR,
    M_dyn_over_Mb_LambdaCDM,
    M_dyn_over_Mb_framework,
)


# ---------------------------------------------------------------------------
# Output paths and matplotlib defaults
# ---------------------------------------------------------------------------

REPO_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
FIG_DIR = os.path.join(REPO_ROOT, "figures")


def _ensure_figdir() -> None:
    os.makedirs(FIG_DIR, exist_ok=True)


def _save(fig, stem: str) -> Tuple[str, str]:
    """Save fig as both PDF and PNG; return (pdf_path, png_path)."""
    _ensure_figdir()
    pdf = os.path.join(FIG_DIR, f"{stem}.pdf")
    png = os.path.join(FIG_DIR, f"{stem}.png")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=200)
    plt.close(fig)
    return pdf, png


def _apply_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "legend.fontsize": 10,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": ":",
    })


# ---------------------------------------------------------------------------
# Figure 1: a_0(z)/a_0(0) = E_eff(z) from z = 0 to z = 15
# ---------------------------------------------------------------------------

def figure1() -> Tuple[str, str, dict]:
    """Figure 1 (§3.2): a_0(z)/a_0(0) vs z, log-uniform y axis.

    Under the anchored convention a_0(z) = a_0_SPARC * E(z)/E(0), this is
    identical to E_eff(z). The paper recommends a log-uniform vertical
    axis; we comply.
    """
    z = np.linspace(0.0, 15.0, 601)
    y = np.array([_E_eff(zi) for zi in z])

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(z, y, color="C0", lw=2.0, label=r"$a_0(z)/a_0(0)\;=\;E(z)/E(0)$")

    # Annotated reference points from paper Table 1.
    ref_z = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0]
    ref_y = [_E_eff(zi) for zi in ref_z]
    ax.scatter(ref_z, ref_y, s=22, color="C3", zorder=5,
               label="Table 1 reference epochs")

    ax.set_yscale("log")
    ax.set_xlabel(r"redshift $z$")
    ax.set_ylabel(r"$a_0(z) / a_0(0)$  [$\equiv E(z)/E(0)$]")
    ax.set_title(r"Figure 1: $a_0(z)/a_0(0)$ from $z=0$ to $z=15$ (flat $\Lambda$CDM, Planck 2018)")
    ax.set_xlim(0.0, 15.0)
    ax.legend(loc="lower right", frameon=False)

    pdf, png = _save(fig, "figure1")

    # ---- Assertion: paper Table 1 has E(2) = 3.033 (3 dp) ----
    val = float(_E_eff(2.0))
    target = 3.033
    assert round(val, 3) == target, (
        f"Figure 1 assertion failed: a_0(2)/a_0(0) = {val:.6f}, "
        f"rounded to 3 dp = {round(val,3)}, paper Table 1 quotes {target}."
    )
    info = {
        "n_points": int(z.size),
        "y_at_z2": val,
        "y_at_z2_rounded": round(val, 3),
        "paper_target_at_z2": target,
        "assertion": "passed",
    }
    return pdf, png, info


# ---------------------------------------------------------------------------
# Figure 2: BTFR loci at z = 0, 1, 2
# ---------------------------------------------------------------------------

def figure2() -> Tuple[str, str, dict]:
    """Figure 2 (§4.3): log10 M_b vs log10 v_flat, three slope-4 lines.

    Anchor: A_BTFR(0) = 62.78 M_sun/(km/s)^4. At redshift z the line is
    log10(M_b) = log10(A_BTFR(0)) - log10(E(z)/E(0)) + 4 log10(v_flat).
    """
    v_kms = np.logspace(np.log10(30.0), np.log10(400.0), 200)
    log_v = np.log10(v_kms)

    A0 = A_BTFR_solar_per_kms4(A0_SPARC_LOCAL)        # 62.78 M_sun/(km/s)^4
    log_A0 = np.log10(A0)

    fig, ax = plt.subplots(figsize=(6.4, 5.0))

    redshifts = [0.0, 1.0, 2.0]
    colors = {0.0: "C0", 1.0: "C2", 2.0: "C3"}
    for z in redshifts:
        Az = A0 * A_ratio(z)                           # = A0 / E_eff(z)
        log_Az = np.log10(Az)
        log_Mb = log_Az + 4.0 * log_v
        ax.plot(log_v, log_Mb,
                color=colors[z], lw=2.0,
                label=fr"$z={z:.0f}$:  $\log_{{10}}A_\mathrm{{BTFR}}={log_Az:.3f}$")

    # Mark archetype points (M_b, v_flat(0)) on the z=0 line as visual anchor.
    for arch in ARCHETYPES:
        v0 = v_flat(arch.M_b, 0.0)                    # km/s
        ax.scatter([np.log10(v0)], [np.log10(arch.M_b)],
                   color="C0", s=28, edgecolor="black",
                   linewidth=0.6, zorder=5)

    ax.set_xlabel(r"$\log_{10}\,v_\text{flat}$ [km/s]")
    ax.set_ylabel(r"$\log_{10}\,M_b$ [$M_\odot$]")
    ax.set_title("Figure 2: BTFR loci at $z=0,1,2$ (slope 4, deep-MOND limit)")
    ax.set_xlim(np.log10(30.0), np.log10(400.0))
    ax.legend(loc="lower right", frameon=False)

    pdf, png = _save(fig, "figure2")

    # ---- Assertion: at log10(v) = log10(176), z=2 line gives log10(M_b)
    # equal to log10(A_BTFR(0)/E(2)) + 4 log10(176). Compare against the
    # §4.2 worked example (M_b = 6e10, v(2) ~= 232 km/s) by checking that
    # the z=2 line evaluated at v = 232 km/s gives log10(M_b) very close
    # to log10(6e10) = 10.7782. ----
    A2 = A0 * A_ratio(2.0)
    log_Mb_at_v2 = np.log10(A2) + 4.0 * np.log10(232.0)
    target = np.log10(6.0e10)
    diff = abs(log_Mb_at_v2 - target)
    assert diff < 0.01, (
        f"Figure 2 assertion failed: z=2 line at v=232 km/s gives "
        f"log10(M_b) = {log_Mb_at_v2:.4f}, paper L* says {target:.4f}; "
        f"|diff| = {diff:.4f} > 0.01 dex."
    )
    info = {
        "v_range_kms": (30.0, 400.0),
        "log10_A_BTFR0": float(log_A0),
        "log10_A_BTFR_z2": float(np.log10(A2)),
        "z=2_line_at_v232_kms": float(log_Mb_at_v2),
        "paper_log10_Mb_Lstar": float(target),
        "abs_diff_dex": float(diff),
        "assertion": "passed",
    }
    return pdf, png, info


# ---------------------------------------------------------------------------
# Figure 3: L* archetype rotation curve at z = 0 and z = 2
# ---------------------------------------------------------------------------

def _v_circ_simple_mond(r_kpc: np.ndarray, M_b_solar: float, z: float) -> np.ndarray:
    """Circular velocity in km/s under the simple MOND interpolation
    mu(x) = x/(1+x), with point-mass Newtonian source g_N(r) = G M_b / r^2.

    Closed form (paper Figure 3 caption):
        g_tot = ( g_N + sqrt(g_N^2 + 4 g_N a_0(z)) ) / 2
        v(r)  = sqrt( r * g_tot )

    All inputs/outputs are in SI internally; r_kpc and the returned
    velocity are in kpc and km/s respectively.
    """
    r_m = r_kpc * M_PER_KPC                                   # m
    M_kg = M_b_solar * M_SUN_KG                               # kg
    a0_z = a0_of_z_anchored(z)                                # m/s^2
    g_N = G_NEWTON * M_kg / r_m ** 2                          # m/s^2
    g_tot = 0.5 * (g_N + np.sqrt(g_N ** 2 + 4.0 * g_N * a0_z))
    v_si = np.sqrt(r_m * g_tot)                               # m/s
    return v_si * KM_PER_M                                    # km/s


def figure3() -> Tuple[str, str, dict]:
    """Figure 3 (§5.3): L* archetype rotation curve at z = 0 and z = 2.

    Per the paper caption: simple MOND mu(x) = x/(1+x) with closed form
    g_tot = (g_N + sqrt(g_N^2 + 4 g_N a_0(z)))/2 and point-mass
    g_N = G M_b / r^2. Plotted only for r >= R_d = 2.5 kpc.
    """
    Lstar = next(a for a in ARCHETYPES if a.name.startswith("L*"))
    R_d = Lstar.R_d                                            # 2.5 kpc

    # Plot to 100 kpc so the slow approach of the simple-MOND v(r) to its
    # BTFR plateau (from above, monotonically) is clearly visible. The
    # caption restricts the meaningful regime to r >= R_d; the upper end
    # is purely a display choice.
    r = np.linspace(R_d, 100.0, 400)                           # kpc

    v_z0 = _v_circ_simple_mond(r, Lstar.M_b, 0.0)
    v_z2 = _v_circ_simple_mond(r, Lstar.M_b, 2.0)

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.plot(r, v_z0, color="C0", lw=2.0, label=r"$z=0$")
    ax.plot(r, v_z2, color="C3", lw=2.0, label=r"$z=2$")

    # Asymptotic plateaus (BTFR) overlay.
    v0_inf = v_flat(Lstar.M_b, 0.0)
    v2_inf = v_flat(Lstar.M_b, 2.0)
    ax.axhline(v0_inf, color="C0", ls=":", alpha=0.6,
               label=fr"$v_\text{{flat}}(0) = {v0_inf:.0f}$ km/s")
    ax.axhline(v2_inf, color="C3", ls=":", alpha=0.6,
               label=fr"$v_\text{{flat}}(2) = {v2_inf:.0f}$ km/s")

    # Mark r_M at each redshift on the curve.
    rM0 = r_M(Lstar.M_b, 0.0)
    rM2 = r_M(Lstar.M_b, 2.0)
    ax.axvline(rM0, color="C0", ls="--", alpha=0.4,
               label=fr"$r_M(0) = {rM0:.2f}$ kpc")
    ax.axvline(rM2, color="C3", ls="--", alpha=0.4,
               label=fr"$r_M(2) = {rM2:.2f}$ kpc")

    ax.set_xlabel(r"radius $r$ [kpc]  ($r \geq R_d = 2.5$ kpc)")
    ax.set_ylabel(r"$v_\mathrm{circ}(r)$ [km/s]")
    ax.set_title(r"Figure 3: $L^*$ rotation curve at $z=0$ and $z=2$ (simple MOND, point-mass)")
    ax.set_xlim(0.0, 100.0)
    ax.set_ylim(150.0, max(v2_inf * 1.18, 270.0))
    ax.legend(loc="upper right", frameon=False, ncol=1)

    pdf, png = _save(fig, "figure3")

    # ---- Assertions ----
    # (i) Simple MOND with point-mass approaches v_flat *from above*
    # monotonically. We check that (a) the outermost plotted v_circ lies
    # above v_flat (small positive overshoot is the physical signature of
    # the slow approach), (b) the curve is monotonically decreasing across
    # the plotted range, (c) the asymptote computed analytically equals
    # the BTFR v_flat, and (d) at the outermost plotted radius the
    # residual Newtonian contribution is below 5%.
    last_v0 = float(v_z0[-1])
    last_v2 = float(v_z2[-1])
    assert last_v0 > v0_inf, (
        f"Figure 3: z=0 outermost v={last_v0:.2f} should sit above the "
        f"v_flat asymptote {v0_inf:.2f} for simple MOND point-mass."
    )
    assert last_v2 > v2_inf, (
        f"Figure 3: z=2 outermost v={last_v2:.2f} should sit above the "
        f"v_flat asymptote {v2_inf:.2f} for simple MOND point-mass."
    )
    assert np.all(np.diff(v_z0) <= 1e-6), (
        "Figure 3: z=0 simple-MOND point-mass v_circ should decrease "
        "monotonically toward v_flat from above."
    )
    assert np.all(np.diff(v_z2) <= 1e-6), (
        "Figure 3: z=2 simple-MOND point-mass v_circ should decrease "
        "monotonically toward v_flat from above."
    )
    rel_excess_0 = (last_v0 - v0_inf) / v0_inf
    rel_excess_2 = (last_v2 - v2_inf) / v2_inf
    assert rel_excess_0 < 0.05, (
        f"Figure 3: residual Newtonian excess at outermost plotted radius "
        f"is {100*rel_excess_0:.2f}% (z=0); should be <5%."
    )
    assert rel_excess_2 < 0.05, (
        f"Figure 3: residual Newtonian excess at outermost plotted radius "
        f"is {100*rel_excess_2:.2f}% (z=2); should be <5%."
    )
    # (ii) Paper plateau values rounded to nearest integer.
    assert round(v0_inf) == 176, (
        f"Figure 3 assertion failed: v_flat(0) = {v0_inf:.4f}, paper 176."
    )
    assert round(v2_inf) == 232, (
        f"Figure 3 assertion failed: v_flat(2) = {v2_inf:.4f}, paper 232."
    )
    info = {
        "M_b_Lstar": Lstar.M_b,
        "R_d_kpc": R_d,
        "r_range_kpc": (float(r[0]), float(r[-1])),
        "v_flat_z0_paper": 176, "v_flat_z0_computed": float(v0_inf),
        "v_flat_z2_paper": 232, "v_flat_z2_computed": float(v2_inf),
        "v_outer_z0": last_v0, "v_outer_z2": last_v2,
        "r_M_z0_kpc": float(rM0), "r_M_z2_kpc": float(rM2),
        "interpolation": "simple MOND mu(x)=x/(1+x); g_tot = (g_N+sqrt(g_N^2+4 g_N a_0))/2",
        "assertion": "passed",
    }
    return pdf, png, info


# ---------------------------------------------------------------------------
# Figure 4: r_M(z) vs z for the four §5 archetypes
# ---------------------------------------------------------------------------

def figure4() -> Tuple[str, str, dict]:
    """Figure 4 (§5.4): r_M(z) for the four §5 archetypes from z = 0 to 5.

    All four curves follow r_M(z) = r_M(0)/sqrt(E(z)/E(0)); they share
    the same fractional shape (mass-independent) and are vertically
    offset by sqrt(G M_b / a_0(0)) for each archetype.
    """
    z = np.linspace(0.0, 5.0, 251)

    fig, ax = plt.subplots(figsize=(6.4, 4.6))

    style = {
        "Dwarf (DDO 154-like)":     ("C0", "-"),
        "Sub-L* (NGC 2403-like)":   ("C2", "-"),
        "L* (NGC 6946-like)":       ("C3", "-"),
        "Giant (UGC 2885-like)":    ("C4", "-"),
    }
    archetype_curves = {}
    for arch in ARCHETYPES:
        rM_curve = np.array([r_M(arch.M_b, zi) for zi in z])
        color, ls = style[arch.name]
        label = (fr"{arch.name}: $r_M(0)={rM_curve[0]:.2f}$ kpc, "
                 fr"$M_b={arch.M_b:.1e}\,M_\odot$")
        ax.plot(z, rM_curve, color=color, ls=ls, lw=2.0, label=label)
        archetype_curves[arch.name] = rM_curve

    ax.set_yscale("log")
    ax.set_xlabel(r"redshift $z$")
    ax.set_ylabel(r"$r_M(z)$ [kpc]")
    ax.set_title(r"Figure 4: MOND transition radius $r_M(z)$ vs $z$ across galaxy mass")
    ax.set_xlim(0.0, 5.0)
    ax.legend(loc="lower left", frameon=False, fontsize=9)

    pdf, png = _save(fig, "figure4")

    # ---- Assertions: §5.2 Table 3 row for L* ----
    Lstar = next(a for a in ARCHETYPES if a.name.startswith("L*"))
    rM0_L = r_M(Lstar.M_b, 0.0)
    rM2_L = r_M(Lstar.M_b, 2.0)
    assert round(rM0_L, 2) == 8.35, (
        f"Figure 4 assertion failed: L* r_M(0) = {rM0_L:.4f}, paper 8.35."
    )
    assert round(rM2_L, 2) == 4.79, (
        f"Figure 4 assertion failed: L* r_M(2) = {rM2_L:.4f}, paper 4.79."
    )
    # Mass-independent fractional shift: r_M(z)/r_M(0) is the same across
    # the archetypes (within float round-off).
    ratios_z2 = [r_M(arch.M_b, 2.0) / r_M(arch.M_b, 0.0) for arch in ARCHETYPES]
    spread = max(ratios_z2) - min(ratios_z2)
    assert spread < 1e-12, (
        f"Figure 4 assertion failed: r_M(2)/r_M(0) is not mass-independent; "
        f"spread = {spread}."
    )
    info = {
        "z_range": (0.0, 5.0),
        "n_archetypes": len(ARCHETYPES),
        "Lstar_rM0_kpc": float(rM0_L),
        "Lstar_rM2_kpc": float(rM2_L),
        "Lstar_rM0_paper": 8.35,
        "Lstar_rM2_paper": 4.79,
        "max_ratio_spread_z2": float(spread),
        "assertion": "passed",
    }
    return pdf, png, info


# ---------------------------------------------------------------------------
# Figure 5: M_dyn(R)/M_b vs R for L* at z = 0, 1, 2 -- framework + ΛCDM
# ---------------------------------------------------------------------------

def figure5() -> Tuple[str, str, dict]:
    """Figure 5 (§6.4): L* M_dyn(R)/M_b vs R for R = 10..500 kpc at
    z = 0, 1, 2.

    Framework prediction: solid lines (eq. 6.3).
    ΛCDM comparison: dashed lines, using Module 5's LSTAR_HALO_TABLE
    (M_halo, c) values per Appendix B.
    """
    R = np.logspace(np.log10(10.0), np.log10(500.0), 200)

    redshifts = (0.0, 1.0, 2.0)
    colors = {0.0: "C0", 1.0: "C2", 2.0: "C3"}

    fig, ax = plt.subplots(figsize=(6.6, 4.8))

    # Framework: solid
    fw_curves = {}
    for z in redshifts:
        y = np.array([M_dyn_over_Mb_framework(Ri, z, M_B_LSTAR) for Ri in R])
        fw_curves[z] = y
        ax.plot(R, y, color=colors[z], ls="-", lw=2.0,
                label=fr"framework $z={z:.0f}$")

    # ΛCDM: dashed (Appendix B halo table)
    lcdm_curves = {}
    for z in redshifts:
        Mh, c = LSTAR_HALO_TABLE[z]
        y = np.array([M_dyn_over_Mb_LambdaCDM(Ri, z, M_B_LSTAR, Mh, c) for Ri in R])
        lcdm_curves[z] = y
        ax.plot(R, y, color=colors[z], ls="--", lw=1.6,
                label=fr"$\Lambda$CDM $z={z:.0f}$ ($M_h={Mh:.1e},\,c={c:.1f}$)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"aperture $R$ [kpc]")
    ax.set_ylabel(r"$M_\mathrm{dyn}(R)\,/\,M_b$")
    ax.set_title(r"Figure 5: $L^*$ ($M_b=6\times10^{10}\,M_\odot$) $M_\text{dyn}/M_b$ vs $R$ at $z=0,1,2$")
    ax.set_xlim(10.0, 500.0)
    ax.legend(loc="upper left", frameon=False, fontsize=9, ncol=2)

    pdf, png = _save(fig, "figure5")

    # ---- Assertions: paper §6.3 quoted values at R = 100 kpc, L* ----
    fw_z0_R100 = M_dyn_over_Mb_framework(100.0, 0.0, M_B_LSTAR)
    fw_z1_R100 = M_dyn_over_Mb_framework(100.0, 1.0, M_B_LSTAR)
    fw_z2_R100 = M_dyn_over_Mb_framework(100.0, 2.0, M_B_LSTAR)
    assert round(fw_z0_R100, 2) == 11.98, (
        f"Figure 5 assertion failed: framework L* z=0 R=100 kpc = "
        f"{fw_z0_R100:.4f}, paper 11.98."
    )
    assert round(fw_z1_R100, 2) == 16.03, (
        f"Figure 5 assertion failed: framework L* z=1 R=100 kpc = "
        f"{fw_z1_R100:.4f}, paper 16.03."
    )
    assert round(fw_z2_R100, 2) == 20.86, (
        f"Figure 5 assertion failed: framework L* z=2 R=100 kpc = "
        f"{fw_z2_R100:.4f}, paper 20.86."
    )
    # ΛCDM at R=100 kpc, paper: 14.0 / 13.6 / 13.8 (1 dp).
    Mh0, c0 = LSTAR_HALO_TABLE[0.0]
    Mh1, c1 = LSTAR_HALO_TABLE[1.0]
    Mh2, c2 = LSTAR_HALO_TABLE[2.0]
    lc_z0_R100 = M_dyn_over_Mb_LambdaCDM(100.0, 0.0, M_B_LSTAR, Mh0, c0)
    lc_z1_R100 = M_dyn_over_Mb_LambdaCDM(100.0, 1.0, M_B_LSTAR, Mh1, c1)
    lc_z2_R100 = M_dyn_over_Mb_LambdaCDM(100.0, 2.0, M_B_LSTAR, Mh2, c2)
    assert round(lc_z0_R100, 1) == 14.0, (
        f"Figure 5 assertion failed: LambdaCDM L* z=0 R=100 kpc = "
        f"{lc_z0_R100:.4f}, paper 14.0."
    )
    assert round(lc_z1_R100, 1) == 13.6, (
        f"Figure 5 assertion failed: LambdaCDM L* z=1 R=100 kpc = "
        f"{lc_z1_R100:.4f}, paper 13.6."
    )
    assert round(lc_z2_R100, 1) == 13.8, (
        f"Figure 5 assertion failed: LambdaCDM L* z=2 R=100 kpc = "
        f"{lc_z2_R100:.4f}, paper 13.8."
    )
    info = {
        "R_range_kpc": (10.0, 500.0),
        "framework_R100_z0": float(fw_z0_R100),
        "framework_R100_z1": float(fw_z1_R100),
        "framework_R100_z2": float(fw_z2_R100),
        "framework_R100_z0_paper": 11.98,
        "lcdm_R100_z0": float(lc_z0_R100),
        "lcdm_R100_z1": float(lc_z1_R100),
        "lcdm_R100_z2": float(lc_z2_R100),
        "assertion": "passed",
    }
    return pdf, png, info


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> None:
    _apply_style()

    print("=" * 70)
    print("Module 11: make_figures.py - generating Figures 1-5")
    print(f"Output directory: {FIG_DIR}")
    print("=" * 70)

    builders = [
        ("Figure 1 (§3.2 a_0(z)/a_0(0))",   figure1),
        ("Figure 2 (§4.3 BTFR loci)",       figure2),
        ("Figure 3 (§5.3 L* rotation curve)", figure3),
        ("Figure 4 (§5.4 r_M vs z)",        figure4),
        ("Figure 5 (§6.4 M_dyn/M_b vs R)",  figure5),
    ]

    results = []
    for name, fn in builders:
        print()
        print(f">>> {name}")
        pdf, png, info = fn()
        print(f"    -> wrote {pdf}")
        print(f"    -> wrote {png}")
        for k, v in info.items():
            print(f"       {k}: {v}")
        results.append((name, pdf, png, info))

    print()
    print("=" * 70)
    print("All 5 figures generated; all assertion checks passed.")
    print("=" * 70)


if __name__ == "__main__":
    main()
