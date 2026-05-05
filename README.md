# a₀(z): Epoch-Dependent MOND Acceleration Scale from Bounded Topology

Code and data for the paper:

**B. Shatto, "Epoch-Dependent Acceleration Scale from Bounded Topology: Predictions for High-Redshift Galactic Dynamics" (2026).**

This repository contains the analysis pipeline, prediction tables, and figure-generation scripts needed to reproduce all numerical results in the paper.

---

## Overview

The paper derives an epoch-dependent MOND acceleration scale

$$a_0(z) = a_0(0) \cdot E(z)$$

within a bounded-topology measurement framework, where $E(z) = H(z)/H_0$ is the Friedmann dimensionless Hubble parameter. The local Milgrom ratio $a_0/(cH_0)$ is predicted as the structural ratio of two phase-operator values at Fibonacci wells of the framework's 120-domain, $C(13/120)/C(34/120) = 0.1845$, agreeing with observation at the 0.7% level.

Five distinct observable channels follow from the single relation $a_0(z) \propto E(z)$ as different powers of $E(z)$:

- BTFR normalization shifts as $E(z)^{-1}$ (§3.1)
- MOND transition radius contracts as $E(z)^{-1/2}$ (§3.2)
- Asymptotic flat velocity at fixed baryonic mass rises as $E(z)^{+1/4}$ (§3.2)
- Lensing-inferred $M_\text{dyn}/M_b$ enhancement scales as $E(z)^{+1/2}$ (§3.3)
- Free-fall collapse-time speedup as $E(z)^{-1/4}$ (§3.4, supporting evidence)

Reproducible results in this repository:

- **a₀(z) tabulation** across $z = 0$ to $15$ (Table 1, §3)
- **BTFR evolution** at six reference redshifts (Table 2, §3.1)
- **Rotation-curve archetype predictions** for four baryonic-mass anchors (Table 3, §3.2)
- **Lensing $M_\text{dyn}/M_b$** at three apertures across five redshifts (Table 4, §3.3)
- **Five-prediction summary** at $z = 2$ (Table 5, §3.5)
- **Constraint-status summary** across five regimes (Table 6, §4.5)
- **Falsification criteria** for each predicted exponent (Table 7, §5)
- **Near-term test schedule** (Table 8, §5)
- **ΛCDM lensing discriminator with pess/repr/opt bracketing** (Table B.1, App B.3)
- **Übler 2017 KMOS3D σ-tension** under three uncertainty budgets (§4.1)
- **Übler 2017 forward-model bias sweep** (four-bias Monte Carlo: RPU, BS, NA, SI) showing no literature-plausible combination reproduces the observed BTFR pattern (§4.1)
- **CMB leakage-coupling bound** $\varepsilon \leq 1.2 \times 10^{-5}$ at 0.5% Planck tolerance (§4.3)
- **Combinatorial baseline** for the Milgrom-ratio match: 24 of 7,021 pairs match within 1%; 6 unique value-pairs after reflection collapse; 1 of 6 within the Fibonacci subset $\{13, 21, 34, 55\}/120$ (§2)
- **Labbé candidate speedup factors** $E(z)^{1/4} \in [1.92, 2.06]$ (§3.4)

---

## Citation

If you use this code or data, please cite the paper and the Zenodo archive of this repository.

The repository is archived on Zenodo:

- **Concept DOI** (resolves to the latest version): [10.5281/zenodo.19980665](https://doi.org/10.5281/zenodo.19980665)
- **Version DOI** (this release, v1.2.0): [10.5281/zenodo.20045115](https://doi.org/10.5281/zenodo.20045115)

For paper citations and machine-readable use the **concept DOI**; it auto-resolves to the latest version. For citing the exact code state used in a specific analysis, use the version DOI.

```bibtex
@software{ShattoA0ZCode2026,
  author = {Shatto, B.},
  title  = {a0z: Code and Data for "Epoch-Dependent Acceleration
            Scale from Bounded Topology"},
  year   = {2026},
  doi    = {10.5281/zenodo.19980665},
  url    = {https://doi.org/10.5281/zenodo.19980665},
  note   = {Concept DOI; resolves to the latest version}
}
```

---

## Repository contents

```
.
├── README.md                            This file
├── LICENSE                              MIT
├── requirements.txt                     Python dependencies
├── data/
│   ├── README.md                        Provenance and reference for observational inputs
│   └── labbe_2023_candidates.csv        Labbé et al. 2023 candidate IDs and photometric redshifts
├── scripts/
│   ├── README.md                        Module overview and verification target list
│   ├── cosmology.py                     E(z), H(z), t_age(z) under flat ΛCDM (Planck 2018)
│   ├── framework.py                     C(Θ), N_H, scaling-law evaluation, well assignments
│   ├── combinatorial_baseline.py        §2 sparsity check on the 120-domain
│   ├── btfr.py                          §3.1 BTFR normalization and velocity scaling
│   ├── rotation_curves.py               §3.2 r_M, v_flat at archetype masses
│   ├── lensing.py                       §3.3 + App B M_dyn/M_b, ΛCDM NFW comparison
│   ├── jwst_speedup.py                  §3.4 free-fall collapse, Labbé candidates
│   ├── ubler_sigma_tension.py           §4.1 analytic σ-tension table
│   ├── ubler_forward_model.py           §4.1 four-bias Monte Carlo forward analysis
│   ├── cmb_leakage.py                   §4.3 ε bound from Planck first-peak amplitude
│   ├── verify.py                        Single-command pass/fail harness across all modules
│   ├── make_tables.py                   Generate Tables 1–8 + Table B.1 from the analysis modules
│   └── make_figures.py                  Generate the four paper figures from the analysis modules
├── figures/
│   ├── figure1.{pdf,png}                §2 Milgrom-ratio sparsity on the 120-domain
│   ├── figure2.{pdf,png}                §3.3 + App B Euclid DR1 lensing discriminator
│   ├── figure3.{pdf,png}                §3.5 five exponents from one input
│   └── figure4.{pdf,png}                §4.1 Übler trend-shape tension
└── tables/
    ├── table{1..8}.csv                  One CSV per paper table (§3 / §3.1 / §3.2 / §3.3 / §3.5 / §4.5 / §5 / §5)
    └── tableB1.csv                      App B.3 ΛCDM bracketing at L*, R = 100 kpc
```

---

## Installation

Python 3.10 or later recommended.

```bash
git clone https://github.com/dmobius3/a0z.git
cd a0z
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Dependencies (`requirements.txt`):

```
numpy>=1.24
scipy>=1.10
matplotlib>=3.7
pandas>=2.0
```

Total install footprint about 100 MB. The full pipeline (all tables + all figures) runs in under a minute on a recent laptop. No MCMC; the predictions are deterministic given the framework's well assignments and the standard ΛCDM cosmology.

---

## Verification

A single command checks every numerical claim in the paper against the corresponding script output:

```bash
python scripts/verify.py
```

Exit code 0 plus `All modules report green against the paper.` means every assertion the paper depends on is reproducible from the scripts. Exit code 1 prints the failure list module-by-module with the relevant stdout tail.

---

## Quickstart: reproducing all paper numbers

```bash
# Generate Tables 1–8 (and Table B.1)
python scripts/make_tables.py

# Generate paper figures
python scripts/make_figures.py
```

Each analysis module also runs standalone and prints its own per-claim verification. Examples:

```bash
python scripts/cosmology.py            # §3 cosmology + Table 1 row spot-checks
python scripts/framework.py            # §2 phase-operator ratios
python scripts/combinatorial_baseline.py  # §2 sparsity (accepts --target, --tol, --list-matches)
python scripts/btfr.py                 # §3.1 Table 2 + L* worked example
python scripts/rotation_curves.py      # §3.2 archetype rotation curves
python scripts/lensing.py              # §3.3 + App B M_dyn/M_b
python scripts/jwst_speedup.py         # §3.4 Labbé speedup factors
python scripts/ubler_sigma_tension.py  # §4.1 σ-tension table
python scripts/ubler_forward_model.py  # §4.1 four-bias forward-model sweep
python scripts/cmb_leakage.py          # §4.3 ε bound
```

Pipeline runtime is well under one minute on a recent laptop. No fitting, no MCMC; predictions are deterministic given the Planck 2018 cosmology and the framework's $z = 0$ well assignments ($\Theta_{a_0} = 13/120$, $\Theta_H = 34/120$).

---

## Data sources and provenance

| Dataset | Source | Reference |
|---|---|---|
| SPARC rotation curves | [astroweb.cwru.edu/SPARC](http://astroweb.cwru.edu/SPARC/) | Lelli, McGaugh, Schombert, *Astron. J.* **152**, 157 (2016) |
| Übler 2017 KMOS3D BTFR | Tabulated from publication | Übler et al., *Astrophys. J.* **842**, 121 (2017) |
| Labbé 2023 candidate sample | [Public Nature catalog](https://www.nature.com/articles/s41586-023-05786-2) | Labbé et al., *Nature* **616**, 266 (2023) |
| Wisnioski 2015 σ_0 medians | KMOS3D survey paper | Wisnioski et al., *Astrophys. J.* **799**, 209 (2015) |
| Planck 2018 cosmology | [pla.esac.esa.int](https://pla.esac.esa.int/) | Planck Collaboration VI, *Astron. Astrophys.* **641**, A6 (2020) |

The Labbé candidate file under `data/` is a formatted derivative of the public Nature supplement above, repackaged for direct loading by `jwst_speedup.py`. No proprietary data is redistributed; the SPARC, Übler, Wisnioski, and Planck values used elsewhere in the pipeline are taken directly from the cited publications and hardcoded in the relevant scripts.

---

## Reproducibility notes

- **No fitting**: the predictions in this paper are deterministic given the framework's well assignments at $z = 0$ ($\Theta_{a_0} = 13/120$, $\Theta_H = 34/120$) and the standard ΛCDM cosmology (Planck 2018). There are no MCMC chains to converge.
- **No Monte Carlo**: every prediction in this pipeline is closed-form. The §4.1 σ-tension table is analytic; the §3.3 / Appendix B ΛCDM comparison uses closed-form NFW enclosed mass; the §3.4 speedup factors are direct evaluations of $E(z)^{1/4}$.
- **Numerical accuracy**: the cosmic-age integral uses SciPy's `quad` with `epsrel=1e-8`; NFW enclosed-mass evaluation is closed-form.
- **Random seeds**: where Monte Carlo is added in subsequent modules (e.g., a forward-model variant of the §4.1 Übler analysis), every script accepts a `--seed` argument; default is 42 for deterministic outputs.
- **Platform**: tested on macOS 14.x and Ubuntu 22.04 with Python 3.11. No GPU required.

---

## License

This repository is released under the MIT License. See `LICENSE` for the full text.

---

## Contact

- Author: B. Shatto, bshatto.pe@gmail.com
- Issues and questions: please open a [GitHub issue](https://github.com/dmobius3/a0z/issues).
