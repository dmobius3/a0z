# Analysis scripts

Python pipeline that reproduces every numerical claim in the paper. To be populated by the script-building unit per the specs in `../README.md`.

## Module overview

| Module | Purpose | Paper section |
|---|---|---|
| `cosmology.py` | Flat ΛCDM E(z), H(z), t_age(z) under Planck 2018 parameters | §3 |
| `framework.py` | C(Θ) phase operator, N_H normalization, well assignments, scaling-law evaluation | §2 |
| `btfr.py` | BTFR normalization shift, Form A/B velocity and mass predictions | §4 |
| `rotation_curves.py` | r_M(z), v_flat(z) at archetype baryonic masses | §5 |
| `lensing.py` | M_dyn/M_b under Newtonian inversion; ΛCDM NFW + Moster + Duffy comparison | §6, Appendix B |
| `jwst_speedup.py` | E(z)^(1/4) free-fall collapse speedup, Labbé candidates | §7.7 |
| `ubler_forward_model.py` | Four-bias-model Monte Carlo simulation of Übler pipeline | §7.2 |
| `cmb_leakage.py` | ε bound from Planck first-peak amplitude precision | §7.5 |
| `combinatorial_baseline.py` | Sparsity check on 120-domain phase-position pairs | §2.5 |
| `make_tables.py` | Generate Tables 1–7 from the analysis modules | All |
| `make_figures.py` | Generate Figures 1–5 from the analysis modules | §§3.2, 4.3, 5.3, 5.4, 6.4 |

## Coding standards

- Python 3.10+
- `numpy`, `scipy`, `matplotlib`, `pandas` only (see `../requirements.txt`)
- `--seed` argument default 42 for any Monte Carlo
- Exact numerical reproducibility against the paper's Tables 1–7
- Closed-form NFW enclosed mass; no halo-mass-function fitting
- Cosmic-age integral via SciPy `quad` with `epsrel=1e-8`

## Verification target numbers

Scripts must hit these values from the paper exactly (within reported precision):

- Phase-operator ratio: $C(13/120)/C(34/120) = 0.1845$
- Observed Milgrom ratio: $a_0/(cH_0) = 0.1833$
- $E(2) = 3.0327$
- BTFR shift at $z = 2$: $1/E(z) = 0.330$
- Lensing enhancement at $z = 2$: $\sqrt{E(z)} = 1.741$
- $L^*$ archetype $r_M(0) = 8.35$ kpc, $v_\text{flat}(0) = 176$ km/s
- ΛCDM L* M_dyn/M_b at R = 100 kpc: $14.0, 13.6, 13.7$ at $z = 0, 1, 2$
- Übler forward model: closest combined model recovers $(-0.295, -0.472)$ dex at $(z=0.9, z=2.3)$
- Σ-tension table: joint $7.6\sigma$ / $5.1\sigma$ / $2.9\sigma$ across three uncertainty budgets
- CMB leakage bound: $\varepsilon \lesssim 1.2 \times 10^{-5}$ at 0.5% Planck tolerance
- Combinatorial sparsity: 24 of 7,021 phase-position pairs match within 1%; 1 of 10 within Fibonacci-well subset
