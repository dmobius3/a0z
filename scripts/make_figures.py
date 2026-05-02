"""make_figures.py
==============

Generate the four figures of the paper from the public APIs of the
analysis modules. Each figure is anchored to a specific paper claim
that benefits from visualization beyond the corresponding table.

Figures produced
----------------

Figure 1 (§2). Milgrom-ratio sparsity on the 120-domain.
    Histogram of all 7,021 unordered phase-operator ratios
    min(C(k1/120)/C(k2/120), C(k2/120)/C(k1/120)) with the 24 matches
    within 1% of 0.1833 highlighted as ticks above the histogram and
    the framework's (13, 34) pair starred. Lands the §2 sparsity
    claim viscerally; the 0.34% number in prose becomes a thicket.

Figure 2 (§3.3 + App B). Euclid DR1 lensing discriminator.
    L* M_dyn(R)/M_b vs aperture R for R = 10..500 kpc at z = 0, 1, 2.
    Framework: solid (universal sqrt(E(z)) shift). ΛCDM: dashed (NFW +
    Moster SHMR + Duffy concentration, App B halo table). At z = 2 the
    ΛCDM curve is bracketed with the App B.1 pessimistic/optimistic
    parameterizations as a shaded band. Lands the §3.3 / App B / §5
    punchline.

Figure 3 (§3.5). Five exponents from one input.
    Five normalized observables vs z on a log-y axis, all driven by
    E(z) under the anchored convention: BTFR normalization (E^-1),
    MOND radius (E^-1/2), v_flat at fixed M_b (E^+1/4), lensing
    M_dyn/M_b (E^+1/2), free-fall collapse time (E^-1/4). Vertical
    dashed line at z = 2 with annotated values matching Table 5
    (0.330, 0.574, 1.320, 1.741, 0.758). The "five exponents from
    one input" structural claim becomes a single image.

Figure 4 (§4.1). Übler trend-shape tension.
    Δb (BTFR zero-point shift) vs z, z = 0..3. Framework prediction:
    clean solid curve = -log10(E(z)/E(0)) under the anchored
    convention; zero free parameters. Übler 2017 KMOS3D points:
    z = 0.9 (Δb_obs = -0.443 ± 0.04), z = 2.3 (-0.27 ± 0.05). Joint
    tension 7.6 / 5.1 / 2.9 σ across the three uncertainty budgets is
    annotated. The asymmetry (clean framework line vs error-barred
    data) is the visual statement: framework has nowhere to hide.

Outputs
-------

Each figure is saved to ../figures/figureN.{pdf,png}.

Modules consumed
----------------

    cosmology              for E(z) and H(z)
    framework              for the phase operator C(theta) and wells
    combinatorial_baseline for the §2 ratio sweep and the 24 matches
    btfr                   for _E_eff(z) (anchored E(z)/E(0))
    rotation_curves        for archetype baryonic masses (informational)
    lensing                for framework + ΛCDM M_dyn/M_b at L*, R
    ubler_sigma_tension    for the Übler 2017 binned points

Each figure carries an in-script self-check against a paper-quoted
value to catch regressions (e.g., the framework value at z = 2,
R = 100 kpc rounds to 20.86; Δb_MIT(0.9) rounds to -0.227).
"""

from __future__ import annotations

import os
import sys
from typing import Iterable, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Ensure sibling modules import whether run from repo root or scripts/.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from cosmology import PLANCK18                                # noqa: E402
from framework import C_phase, DOMAIN_SIZE                    # noqa: E402
from combinatorial_baseline import (                          # noqa: E402
    OBSERVED_RATIO,
    TOLERANCE_FRAC,
    smaller_ratio,
    matches_within,
    canonical_under_reflection,
)
from btfr import _E_eff                                       # noqa: E402
from lensing import (                                         # noqa: E402
    M_B_LSTAR,
    M_dyn_over_Mb_framework,
    M_dyn_over_Mb_LambdaCDM,
    LSTAR_HALO_TABLE,
)
from ubler_sigma_tension import UBLER_BINS, delta_b_MIT       # noqa: E402


# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------

FIGURES_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "figures"))


def _ensure_figures_dir() -> None:
    os.makedirs(FIGURES_DIR, exist_ok=True)


def _save(fig: plt.Figure, stem: str) -> Tuple[str, str]:
    """Save fig to figures/<stem>.pdf and figures/<stem>.png. Return paths."""
    pdf = os.path.join(FIGURES_DIR, f"{stem}.pdf")
    png = os.path.join(FIGURES_DIR, f"{stem}.png")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return pdf, png


# ---------------------------------------------------------------------------
# Figure 1 (§2): Milgrom-ratio sparsity on the 120-domain
# ---------------------------------------------------------------------------

def figure1_sparsity() -> Tuple[str, str]:
    """Figure 1 (§2): two-panel sparsity figure on the 120-domain.

    NOTE: This auto-generator is kept for reference and regeneration but
    is NOT called from main(). The shipped figures/figure1.{pdf,png} is
    a curated artifact (the (13, 34) tag was removed for layout cleanup
    while preserving the legend label and the rest of the figure
    content). Calling this function will overwrite the curated PDF/PNG
    with the auto version. Do that only when intentionally regenerating
    from scratch.

    Top: histogram of all 7,021 unordered phase-operator ratios in [0, 1],
    with a vertical line at the observed Milgrom ratio 0.1833 and a
    shaded ±1% band marking the matching region.

    Bottom: zoom into the band region [0.176, 0.190] showing each of the
    24 matching ratios as a distinct tick at its actual x-position, with
    the framework's (13, 34) Fibonacci pair starred and labeled. The
    bottom panel makes the 24 sparse matches individually visible.
    """
    target = OBSERVED_RATIO
    tol = TOLERANCE_FRAC
    band_lo = target * (1 - tol)
    band_hi = target * (1 + tol)

    # All 7,021 unordered nonzero ratios on the 120-domain.
    indices = list(range(1, DOMAIN_SIZE))
    all_ratios: List[float] = []
    matching_pairs: List[Tuple[int, int, float]] = []
    for i in range(len(indices)):
        for j in range(i + 1, len(indices)):
            k1, k2 = indices[i], indices[j]
            r = smaller_ratio(k1, k2)
            all_ratios.append(r)
            if matches_within(k1, k2, target=target, tol=tol):
                matching_pairs.append((k1, k2, r))

    fig, (ax_top, ax_bot) = plt.subplots(
        nrows=2, ncols=1, figsize=(7.4, 5.6),
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.42},
    )

    # ---- top panel: full distribution ----
    bins = np.linspace(0.0, 1.0, 80)
    ax_top.hist(
        all_ratios, bins=bins, color="lightgrey",
        edgecolor="grey", linewidth=0.4, alpha=0.85,
        label=r"All 7,021 unordered pairs on $S^3/2I$",
    )
    ax_top.axvspan(
        band_lo, band_hi, color="tab:red", alpha=0.30,
        label=r"$\pm 1\%$ of observed $a_0/(cH_0) = 0.1833$",
    )
    ax_top.axvline(target, color="tab:red", linestyle="--", linewidth=1.2)
    ax_top.set_xlabel(
        r"$\min\,[\,C(k_1/120)/C(k_2/120),\ C(k_2/120)/C(k_1/120)\,]$"
    )
    ax_top.set_ylabel("Number of unordered pairs")
    ax_top.set_title(
        r"Figure 1 (§2): Milgrom-ratio sparsity on the 120-domain"
    )
    ax_top.set_xlim(0.0, 1.0)
    ax_top.legend(loc="upper right", frameon=False, fontsize=9)

    # ---- bottom panel: zoom into the band, ticks for each match ----
    zoom_lo, zoom_hi = 0.176, 0.190
    ax_bot.axvspan(band_lo, band_hi, color="tab:red", alpha=0.18)
    ax_bot.axvline(target, color="tab:red", linestyle="--", linewidth=1.2)
    ax_bot.text(
        target, 1.07, fr"$0.1833$",
        color="tab:red", fontsize=9, ha="center", va="bottom",
    )

    # Each match as a vertical tick at its actual ratio. Star (13, 34).
    for k1, k2, r in matching_pairs:
        is_framework = (k1, k2) == (13, 34) or (k1, k2) == (34, 13)
        if is_framework:
            ax_bot.scatter(
                [r], [0.50], marker="*", s=180,
                color="tab:blue", edgecolor="black", linewidth=0.6,
                zorder=5,
            )
            # Tag the (13, 34) star above the right-side caption box so
            # that both annotations live on the right-hand side. The
            # leader is an L-shape (angle3 connection): leaves the text
            # going LEFT, then turns and approaches the star from
            # ABOVE. The horizontal leg sits at y ≈ 1.05 (above box
            # top); the vertical leg is at x = r (left of the box).
            # The line never crosses the caption-box interior.
            ax_bot.annotate(
                r"$(13,\,34)$", xy=(r, 0.50),
                xytext=(zoom_hi - 0.0005, 1.05),
                fontsize=10, color="tab:blue",
                ha="right", va="center",
                arrowprops=dict(
                    arrowstyle="-", color="tab:blue", linewidth=0.6,
                    connectionstyle="angle3,angleA=90,angleB=180",
                ),
            )
        else:
            ax_bot.vlines(r, 0.20, 0.80, color="tab:orange",
                          linewidth=1.6, zorder=3)

    # Legend strokes for the bottom panel.
    ax_bot.vlines([], [], [], color="tab:orange", linewidth=1.6,
                  label=fr"{len(matching_pairs)} matching pairs (within $1\%$)")
    ax_bot.scatter([], [], marker="*", s=180, color="tab:blue",
                   edgecolor="black", linewidth=0.6,
                   label=r"Framework Fibonacci pair $(13,\,34)$")
    ax_bot.legend(loc="lower left", frameon=False, fontsize=8.5)

    # Cosmetics for the bottom panel.
    ax_bot.set_xlim(zoom_lo, zoom_hi)
    ax_bot.set_ylim(0.0, 1.15)
    ax_bot.set_xlabel("Phase-operator ratio (zoom into the matching band)")
    ax_bot.set_yticks([])
    for spine_name in ("left", "top", "right"):
        ax_bot.spines[spine_name].set_visible(False)

    # Caption-style note on the right of the bottom panel.
    ax_bot.text(
        zoom_hi - 0.0005, 0.5,
        f"  24 of 7,021 = 0.34%\n"
        f"  6 unique value-pairs\n"
        f"  (mod reflection)\n"
        f"  1 of 6 Fibonacci",
        fontsize=8.5, color="black", ha="right", va="center",
        bbox=dict(boxstyle="round,pad=0.35", fc="white",
                  ec="grey", alpha=0.9),
    )

    # Self-check: count of matches and uniqueness of (13, 34) Fibonacci.
    assert len(matching_pairs) == 24, (
        f"Figure 1: expected 24 matching pairs, got {len(matching_pairs)}"
    )
    canonicals = {canonical_under_reflection(k1, k2) for k1, k2, _ in matching_pairs}
    assert len(canonicals) == 6, (
        f"Figure 1: expected 6 unique value-pairs after reflection, "
        f"got {len(canonicals)}"
    )

    return _save(fig, "figure1")


# ---------------------------------------------------------------------------
# Figure 2 (§3.3 + App B): Euclid DR1 lensing discriminator
# ---------------------------------------------------------------------------

# App B.1 pess/repr/opt halo parameterizations at L*, used for the bracketing
# band at z = 2 (and z = 0 if we ever want to extend).
APP_B1_HALO = {
    "pessimistic":   {0.0: (2.0e12, 9.0), 2.0: (1.0e12, 4.5)},
    "representative": {0.0: (1.5e12, 7.5), 2.0: (7.0e11, 3.5)},
    "optimistic":    {0.0: (1.0e12, 6.0), 2.0: (5.0e11, 2.5)},
}


def _lcdm_curve(R_kpc: np.ndarray, z: float, M_halo: float, c: float) -> np.ndarray:
    return np.array([
        M_dyn_over_Mb_LambdaCDM(R, z, M_B_LSTAR, M_halo, c) for R in R_kpc
    ])


def _framework_curve(R_kpc: np.ndarray, z: float) -> np.ndarray:
    return np.array([
        M_dyn_over_Mb_framework(R, z, M_B_LSTAR) for R in R_kpc
    ])


def figure2_lensing_discriminator() -> Tuple[str, str]:
    """Figure 2 (§3.3 + App B): L* M_dyn(R)/M_b vs aperture R for
    R = 10..500 kpc at z = 0, 1, 2. Framework solid; ΛCDM dashed
    (representative App B halo table); z = 2 ΛCDM bracketed by
    pess/opt shading from App B.1.
    """
    R = np.linspace(10.0, 500.0, 80)

    fig, ax = plt.subplots(figsize=(7.2, 5.0))

    z_list = [0.0, 1.0, 2.0]
    colors = {0.0: "tab:blue", 1.0: "tab:green", 2.0: "tab:red"}

    # Framework: solid curves (universal sqrt(E(z)) shift).
    for z in z_list:
        ax.plot(R, _framework_curve(R, z),
                color=colors[z], linewidth=2.0,
                label=fr"Framework, $z = {z:g}$")

    # ΛCDM at representative parameters: dashed curves.
    for z in z_list:
        M_h, c = LSTAR_HALO_TABLE[z]
        ax.plot(R, _lcdm_curve(R, z, M_h, c),
                color=colors[z], linewidth=1.6, linestyle="--",
                label=fr"$\Lambda$CDM, $z = {z:g}$ (representative)")

    # Shaded bracketing band at z = 2 between pess and opt ΛCDM curves.
    z_brk = 2.0
    M_h_pess, c_pess = APP_B1_HALO["pessimistic"][z_brk]
    M_h_opt, c_opt = APP_B1_HALO["optimistic"][z_brk]
    pess_curve = _lcdm_curve(R, z_brk, M_h_pess, c_pess)
    opt_curve = _lcdm_curve(R, z_brk, M_h_opt, c_opt)
    band_lo = np.minimum(pess_curve, opt_curve)
    band_hi = np.maximum(pess_curve, opt_curve)
    ax.fill_between(R, band_lo, band_hi,
                    color=colors[z_brk], alpha=0.13,
                    label=r"App B.1 pess/opt bracketing at $z = 2$")

    ax.set_xlabel(r"Aperture $R$ [kpc]")
    ax.set_ylabel(r"$M_\mathrm{dyn}(R)\,/\,M_b$ at $L^*$")
    ax.set_title(r"Figure 2 (§3.3 + App B): Euclid DR1 lensing discriminator")
    ax.set_xlim(10, 500)
    ax.legend(loc="upper left", frameon=False, fontsize=8.5, ncol=1)
    ax.grid(True, which="both", linestyle=":", alpha=0.4)

    # Self-check against paper-quoted values at L*, R = 100 kpc.
    fw_z2_R100 = M_dyn_over_Mb_framework(100.0, 2.0, M_B_LSTAR)
    assert abs(round(fw_z2_R100, 2) - 20.86) < 1e-9, (
        f"Figure 2: framework L* @ R=100, z=2 = {fw_z2_R100:.3f}, "
        f"expected 20.86"
    )
    M_h_repr, c_repr = LSTAR_HALO_TABLE[2.0]
    lcdm_z2_R100 = M_dyn_over_Mb_LambdaCDM(100.0, 2.0, M_B_LSTAR,
                                            M_h_repr, c_repr)
    assert abs(round(lcdm_z2_R100, 1) - 13.8) < 1e-9, (
        f"Figure 2: ΛCDM L* @ R=100, z=2 (repr) = {lcdm_z2_R100:.3f}, "
        f"expected 13.8"
    )

    return _save(fig, "figure2")


# ---------------------------------------------------------------------------
# Figure 3 (§3.5): Five exponents from one input
# ---------------------------------------------------------------------------

# Each row: (label, exponent, color). The observable equals E(z)^exponent
# under the anchored convention E(0) = 1.
FIVE_EXPONENTS = [
    (r"BTFR norm $A(z)/A(0) = 1/E(z)$",       -1.0,    "tab:red"),
    (r"MOND radius $r_M(z)/r_M(0)$",          -0.5,    "tab:orange"),
    (r"$v_\mathrm{flat}$ at fixed $M_b$",      0.25,   "tab:green"),
    (r"Lensing $M_\mathrm{dyn}/M_b$",          0.5,    "tab:blue"),
    (r"Collapse time $t_\mathrm{ff}(z)/t_\mathrm{ff}(0)$",  -0.25, "tab:purple"),
]


def figure3_five_exponents() -> Tuple[str, str]:
    """Figure 3 (§3.5): five normalized observables vs z, all driven by
    a single anchored E(z). Vertical dashed line at z = 2 with annotated
    values matching Table 5 (0.330, 0.574, 1.320, 1.741, 0.758).
    """
    z_grid = np.linspace(0.0, 5.0, 200)
    E_anchored = np.array([_E_eff(z) for z in z_grid])

    fig, ax = plt.subplots(figsize=(7.2, 5.2))

    table5_at_z2: List[float] = []
    for label, exponent, color in FIVE_EXPONENTS:
        y = E_anchored ** exponent
        ax.plot(z_grid, y, color=color, linewidth=1.8, label=label)
        # Tag exponent at right edge.
        ax.text(z_grid[-1] * 1.005, y[-1],
                fr" $E^{{{exponent:+g}}}$",
                color=color, fontsize=9, va="center")
        # z = 2 anchor value (Table 5 entry).
        z2_val = _E_eff(2.0) ** exponent
        table5_at_z2.append(z2_val)

    # Vertical dashed line at z = 2 with annotated values.
    ax.axvline(2.0, color="grey", linestyle=":", alpha=0.7)
    annotation = "  Table 5 ($z = 2$):\n"
    for (label, exp, color), val in zip(FIVE_EXPONENTS, table5_at_z2):
        annotation += f"  $E^{{{exp:+g}}}$ = {val:.3f}\n"
    ax.text(2.05, 0.40, annotation.rstrip(),
            fontsize=8.5, color="black",
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec="grey", alpha=0.9))

    ax.set_yscale("log")
    ax.set_xlabel(r"Redshift $z$")
    ax.set_ylabel(r"Observable / observable at $z = 0$")
    ax.set_title(r"Figure 3 (§3.5): Five exponents from one input")
    ax.set_xlim(0.0, 5.2)
    ax.set_ylim(0.04, 30.0)
    ax.legend(loc="upper left", frameon=False, fontsize=8.5)
    ax.grid(True, which="both", linestyle=":", alpha=0.4)

    # Self-check: round to 3 dp at z = 2 against paper Table 5.
    expected_z2 = [0.330, 0.574, 1.320, 1.741, 0.758]
    rounded = [round(v, 3) for v in table5_at_z2]
    assert rounded == expected_z2, (
        f"Figure 3: Table 5 values at z=2 = {rounded}, "
        f"expected {expected_z2}"
    )

    return _save(fig, "figure3")


# ---------------------------------------------------------------------------
# Figure 4 (§4.1): Übler trend-shape tension
# ---------------------------------------------------------------------------

def figure4_ubler_tension() -> Tuple[str, str]:
    """Figure 4 (§4.1): Δb (BTFR zero-point shift) vs z. Framework
    monotonic curve = -log10(E(z)/E(0)); Übler 2017 KMOS3D points with
    stat-only error bars. Joint tension annotated.
    """
    z_grid = np.linspace(0.0, 3.0, 200)
    framework_curve = np.array([delta_b_MIT(z) for z in z_grid])

    fig, ax = plt.subplots(figsize=(7.0, 4.8))

    # Framework prediction: clean solid line, zero free parameters.
    ax.plot(z_grid, framework_curve,
            color="tab:blue", linewidth=2.0,
            label=r"Framework: $\Delta b_\mathrm{MIT}(z) = -\log_{10}(E(z)/E(0))$")

    # Übler 2017 binned points.
    z_pts = sorted(UBLER_BINS.keys())
    delta_obs = np.array([UBLER_BINS[z][0] for z in z_pts])
    sigma_stat = np.array([UBLER_BINS[z][1] for z in z_pts])
    ax.errorbar(z_pts, delta_obs, yerr=sigma_stat,
                fmt="o", color="tab:red", markersize=8,
                ecolor="tab:red", elinewidth=1.6, capsize=4,
                label=r"Übler+ 2017 KMOS3D (stat only)")

    # Annotate each data point with its Δb_obs value.
    for z, d, s in zip(z_pts, delta_obs, sigma_stat):
        ax.annotate(fr"  $\Delta b = {d:+.3f} \pm {s:.2f}$",
                    xy=(z, d), xytext=(z + 0.08, d - 0.04),
                    fontsize=8.5, color="tab:red")

    # Tension annotation box.
    tension_text = (
        "Joint tension across budgets:\n"
        r"  stat only:        $7.6\sigma$" + "\n"
        r"  $+$ Lelli (0.05 dex):       $5.1\sigma$" + "\n"
        r"  $+$ velcorr (0.10 dex):  $2.9\sigma$"
    )
    ax.text(0.02, 0.04, tension_text,
            transform=ax.transAxes, fontsize=8.5,
            verticalalignment="bottom",
            bbox=dict(boxstyle="round,pad=0.4", fc="white",
                      ec="grey", alpha=0.9))

    # Trend-shape note.
    ax.text(2.5, -0.08, "framework: monotonic\ndata: non-monotonic",
            fontsize=9, color="black", ha="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow",
                      ec="grey", alpha=0.85))

    ax.axhline(0.0, color="grey", linestyle=":", alpha=0.5)
    ax.set_xlabel(r"Redshift $z$")
    ax.set_ylabel(r"BTFR zero-point shift $\Delta b$ [dex]")
    ax.set_title(r"Figure 4 (§4.1): Übler trend-shape tension")
    ax.set_xlim(0.0, 3.0)
    ax.set_ylim(-0.65, 0.10)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.grid(True, linestyle=":", alpha=0.4)

    # Self-check: framework values at z = 0.9, 2.3 round to paper values.
    d_09 = round(delta_b_MIT(0.9), 3)
    d_23 = round(delta_b_MIT(2.3), 3)
    assert d_09 == -0.227, f"Figure 4: Δb_MIT(0.9) = {d_09}, expected -0.227"
    assert d_23 == -0.540, f"Figure 4: Δb_MIT(2.3) = {d_23}, expected -0.540"

    return _save(fig, "figure4")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> int:
    _ensure_figures_dir()

    # Figure 1 (§2 sparsity) is curated externally and shipped as a
    # static PDF/PNG in ../figures/figure1.{pdf,png}. The auto-generator
    # `figure1_sparsity()` is preserved for reference / regeneration but
    # is NOT called from main() so a default `python make_figures.py`
    # run does not clobber the curated artifact. To regenerate the
    # auto version, call `figure1_sparsity()` directly or pass --regen-1.
    #
    # Figures 2, 3, 4 are auto-generated and verified against paper
    # values on every run.
    builders = [
        ("Figure 2 (§3.3 + App B lensing)",       figure2_lensing_discriminator),
        ("Figure 3 (§3.5 five exponents)",        figure3_five_exponents),
        ("Figure 4 (§4.1 Übler tension)",         figure4_ubler_tension),
    ]

    print()
    print("=" * 70)
    print(" Generating paper figures (Figure 1 curated externally; skipped)")
    print("=" * 70)
    print()

    for label, fn in builders:
        try:
            pdf, png = fn()
            print(f"  [OK]  {label}")
            print(f"        -> {pdf}")
            print(f"        -> {png}")
        except AssertionError as e:
            print(f"  [FAIL]  {label}")
            print(f"          self-check failed: {e}")
            return 1
        except Exception as e:
            print(f"  [FAIL]  {label}")
            print(f"          {type(e).__name__}: {e}")
            return 1

    print()
    print("Auto-generated figures (2, 3, 4) pass self-checks. Figure 1 ")
    print("is shipped as a curated PDF/PNG and not regenerated by default.")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
