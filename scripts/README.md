# Analysis scripts

Python pipeline that reproduces every numerical claim in the paper. The lean PRD draft's tables and headline numbers are regenerated cell-by-cell from the modules below; `verify.py` provides a single-command pass/fail gate.

## Module overview

| Module | Purpose | Paper section |
|---|---|---|
| `cosmology.py` | Flat ΛCDM E(z), H(z), t_age(z) under Planck 2018 parameters | §3 (Table 1) |
| `framework.py` | C(Θ) phase operator, N_H normalization, well assignments, scaling-law evaluation | §2 |
| `combinatorial_baseline.py` | Sparsity check on 120-domain phase-position pairs | §2 |
| `btfr.py` | BTFR normalization shift, velocity and mass predictions | §3.1 (Table 2) |
| `rotation_curves.py` | r_M(z), v_flat(z) at archetype baryonic masses | §3.2 (Table 3) |
| `lensing.py` | M_dyn/M_b under Newtonian inversion; ΛCDM NFW + Moster + Duffy comparison | §3.3, App B (Tables 4, B.1) |
| `jwst_speedup.py` | E(z)^(1/4) free-fall collapse speedup, Labbé candidates | §3.4 |
| `ubler_sigma_tension.py` | Analytic σ-tension table for the Übler+ 2017 KMOS3D BTFR comparison | §4.1 |
| `cmb_leakage.py` | ε bound from Planck first-peak amplitude precision | §4.3 |
| `verify.py` | Single-command pass/fail harness across all modules | All |
| `make_tables.py` | Generate Tables 1–8 + Table B.1 from the analysis modules | All |
| `make_figures.py` | Generate paper figures from the analysis modules | TBD |

## Coding standards

- Python 3.10+
- `numpy`, `scipy`, `matplotlib`, `pandas` only (see `../requirements.txt`)
- `--seed` argument default 42 for any Monte Carlo
- Exact numerical reproducibility against the paper's Tables 1–8 and Table B.1
- Closed-form NFW enclosed mass; no halo-mass-function fitting
- Cosmic-age integral via SciPy `quad` with `epsrel=1e-8`

## Verification target numbers

Scripts must hit these values from the paper exactly (within reported precision):

- Phase-operator ratio: $C(13/120)/C(34/120) = 0.1845$
- Observed Milgrom ratio: $a_0/(cH_0) = 0.1833$
- Predicted-vs-observed agreement: $0.7\%$
- $E(2) = 3.033$ (4-dp script value 3.0327)
- $\rho_\text{crit}(0) = 8.53 \times 10^{-27}$ kg/m³
- BTFR shift at $z = 2$: $1/E(z) = 0.330$
- Lensing enhancement at $z = 2$: $\sqrt{E(z)} = 1.741$
- $L^*$ archetype $r_M(0) = 8.35$ kpc, $v_\text{flat}(0) = 176$ km/s
- ΛCDM L* M_dyn/M_b at R = 100 kpc: $14.0, 13.6, 13.8$ at $z = 0, 1, 2$
- App B bracketing (framework / ΛCDM at z=2, R=100 kpc): pessimistic $1.19$ / representative $1.52$ / optimistic $1.86$
- Übler residuals (observed minus framework prediction): $-0.216$ dex at $z=0.9$, $+0.270$ dex at $z=2.3$
- σ-tension table per-bin: $T(z=0.9) = -5.4 / -3.4 / -1.8$ σ; $T(z=2.3) = +5.4 / +3.8 / +2.2$ σ across the three budgets
- σ-tension table joint: $7.6\sigma$ / $5.1\sigma$ / $2.9\sigma$ across three uncertainty budgets
- CMB leakage bound: $\varepsilon \leq 1.2 \times 10^{-5}$ at 0.5% Planck tolerance
- JWST speedup: $E(z)^{1/4} \in [1.92, 2.06]$ across six Labbé candidates
- Combinatorial sparsity: 24 of 7,021 phase-position pairs match within 1%; 6 unique phase-operator value-pairs after the reflection collapse $C(k) = C(120-k)$; 1 of 6 within the Fibonacci-well subset $\{13, 21, 34, 55\}/120$

`verify.py` asserts every value above against the relevant module's API and exits 0 on green. Per-module `main()` blocks also self-verify on direct invocation.
