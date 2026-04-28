# a0z: Epoch-Dependent MOND Acceleration Scale from Bounded Topology

Code and data for the paper:

**B. Shatto, "Epoch-Dependent Acceleration Scale from Bounded Topology: Predictions for High-Redshift Galactic Dynamics" (2026).** Submitted to *Physical Review D*.

This repository contains the analysis pipeline, prediction tables, and figure-generation scripts needed to reproduce all numerical results in the paper.

---

## Overview

The paper derives an epoch-dependent MOND acceleration scale

$$a_0(z) = a_0(0) \cdot E(z)$$

within Mode Identity Theory, where $E(z) = H(z)/H_0$ is the standard Friedmann dimensionless Hubble parameter. The local Milgrom ratio $a_0/(cH_0)$ is predicted as the structural ratio of two phase-operator values at adjacent Fibonacci wells of the framework's 120-domain, $C(13/120)/C(34/120) = 0.1846$, agreeing with observation at the 0.8% level.

Five distinct observable channels follow from the single relation $a_0(z) \propto E(z)$ as different powers of $E(z)$:

- BTFR normalization shifts as $E(z)^{-1}$ (§4)
- MOND transition radius contracts as $E(z)^{-1/2}$ (§5)
- Asymptotic flat velocity at fixed baryonic mass rises as $E(z)^{+1/4}$ (§4–5)
- Lensing-inferred $M_\text{dyn}/M_b$ enhancement scales as $E(z)^{+1/2}$ (§6)
- Free-fall collapse-time speedup as $E(z)^{-1/4}$ (§7.7, supporting evidence)

Reproducible results in this repository:

- **a₀(z) tabulation** across $z = 0$ to $15$ (Table 1, §3.1)
- **BTFR evolution** at six reference redshifts (Table 2, §4.2)
- **Rotation-curve archetype predictions** for four baryonic-mass anchors (Table 3, §5.2)
- **Lensing $M_\text{dyn}/M_b$** at three apertures across five redshifts (Table 4, §6.3)
- **ΛCDM lensing discriminator** (NFW + Moster SHMR + Duffy concentration; §6.2 and Appendix B)
- **Übler 2017 forward-model analysis** across four bias models (§7.2)
- **CMB leakage-coupling bound** $\varepsilon \lesssim 10^{-5}$ from Planck first-peak amplitude (§7.5)
- **Combinatorial baseline** for the Milgrom-ratio match (§2.5)
- **Übler $\sigma$-tension table** under three uncertainty budgets (§7.2)
- **Labbé candidate speedup factors** for the §7.7 supporting analysis
- **Constraint summary** (Table 5), **predictions summary** (Table 6), **falsification criteria** (Table 7)

---

## Citation

If you use this code or data, please cite the paper and the Zenodo archive of this repository.

```bibtex
@article{Shatto2026A0Z,
  title   = {Epoch-Dependent Acceleration Scale from Bounded Topology:
             Predictions for High-Redshift Galactic Dynamics},
  author  = {Shatto, B.},
  journal = {Physical Review D},
  year    = {2026},
  note    = {Submitted}
}

@misc{ShattoA0ZCode2026,
  author = {Shatto, B.},
  title  = {a0z: Code and Data for "Epoch-Dependent Acceleration Scale
            from Bounded Topology"},
  year   = {2026},
  doi    = {[Zenodo DOI to insert at submission]}
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
│   ├── ubler_2017_btfr.csv              Übler et al. 2017 BTFR zero-points and uncertainties
│   ├── labbe_2023_candidates.csv        Labbé et al. 2023 candidate IDs, photo-z, stellar masses
│   └── lelli_2016_sparc_summary.csv     SPARC sample summary (local BTFR baseline)
├── scripts/
│   ├── cosmology.py                     E(z), H(z), t_age(z) under flat ΛCDM (Planck 2018)
│   ├── framework.py                     C(Θ), N_H, scaling-law evaluation, well assignments
│   ├── btfr.py                          §4 BTFR normalization and velocity scaling
│   ├── rotation_curves.py               §5 r_M, v_flat at archetype masses
│   ├── lensing.py                       §6 M_dyn/M_b, NFW comparison (Appendix B)
│   ├── jwst_speedup.py                  §7.7 free-fall collapse, Labbé candidates
│   ├── ubler_forward_model.py           §7.2 four-bias-model forward simulation
│   ├── cmb_leakage.py                   §7.5 ε bound from Planck first-peak amplitude
│   ├── combinatorial_baseline.py        §2.5 sparsity check on 120-domain
│   ├── make_tables.py                   Generate Tables 1–7 from the above
│   └── make_figures.py                  Generate Figures 1–5 from the above
├── figures/
│   ├── fig1_a0_evolution.pdf            §3.2: a_0(z)/a_0(0) vs z
│   ├── fig2_btfr_loci.pdf               §4.3: BTFR loci at z = 0, 1, 2
│   ├── fig3_lstar_rotation.pdf          §5.3: L* rotation curve at z = 0, 2
│   ├── fig4_rm_vs_z.pdf                 §5.4: r_M vs z across four archetypes
│   └── fig5_mdyn_aperture.pdf           §6.4: M_dyn/M_b vs aperture, ΛCDM overlay
├── tables/
│   ├── table1_a0_evolution.csv          §3.1
│   ├── table2_btfr.csv                  §4.2
│   ├── table3_rotation_curves.csv       §5.2
│   ├── table4_lensing.csv               §6.3
│   ├── table5_constraint_summary.csv    §7.8
│   ├── table6_predictions.csv           §9.1
│   └── table7_falsification.csv         §9.2
└── results/                             Supplementary intermediate outputs
    ├── ubler_forward_model.csv          Four bias-model recovered offsets
    ├── lcdm_discriminator_lstar.csv     L* ΛCDM M_dyn/M_b at three apertures × five redshifts
    ├── lcdm_aperture_overlay.csv        60-point aperture grid, ΛCDM curves z = 0, 1, 2
    ├── cmb_leakage_bound.csv            Planck constraint on ε across scales
    ├── combinatorial_baseline.csv       7,021 phase-position pairs and match flags
    ├── ubler_sigma_tension.csv          Three uncertainty budgets (statistical, +Lelli, +velocity)
    └── labbe_candidate_speedup.csv      E(z)^{1/4} per Labbé candidate
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

## Quickstart: reproducing all paper numbers

```bash
# Generate Tables 1–7
python scripts/make_tables.py --output tables/

# Generate Figures 1–5
python scripts/make_figures.py --output figures/
```

Each table CSV and each figure PDF is regenerated from the underlying scripts. Expected median runtimes on a recent laptop: tables under 5 seconds, figures under 30 seconds.

To reproduce a single table or figure:

```bash
# Just the a_0(z) evolution table
python scripts/make_tables.py --table 1 --output tables/table1_a0_evolution.csv

# Just the M_dyn/M_b vs aperture figure
python scripts/make_figures.py --figure 5 --output figures/fig5_mdyn_aperture.pdf
```

---

## Reproducing individual analyses

### §2.5 combinatorial baseline (Milgrom-ratio sparsity)

```bash
python scripts/combinatorial_baseline.py \
    --target-ratio 0.1832 \
    --tolerance 0.01 \
    --output results/combinatorial_baseline.csv
```

Sweeps all 7,021 unordered nonzero phase-position pairs $(k_1, k_2)/120$ with $k_1, k_2 \in \{1, \ldots, 119\}$, computes $C(k_1/120)/C(k_2/120)$, and flags pairs matching the observed Milgrom ratio within 1%. Should report 24 matches on the full 120-domain and 1 unique match within the framework's selected Fibonacci-well subset $\{13, 21, 34, 55, 60\}/120$ (the framework's $(13, 34)$ pair).

### §6.2 ΛCDM lensing discriminator (Appendix B methodology)

```bash
python scripts/lensing.py --discriminator \
    --archetype lstar \
    --apertures 30 100 300 \
    --redshifts 0 0.5 1 2 5 \
    --output results/lcdm_discriminator_lstar.csv
```

Computes both the framework $M_\text{dyn}/M_b = (R/r_M(0)) \sqrt{E(z)}$ and the ΛCDM NFW prediction at each archetype × aperture × redshift combination, using the representative $(M_\text{halo}, c)$ values specified in Appendix B.

### §7.2 Übler forward-model analysis

```bash
python scripts/ubler_forward_model.py \
    --bias-models radial-position beam-smearing non-asymptotic mass-distribution \
    --redshifts 0.9 2.3 \
    --output results/ubler_forward_model.csv
```

Generates mock galaxy samples at $z = 0.9$ and $z = 2.3$ under the framework's BTFR (4.3), applies Übler's radial pressure-support correction with sample-realistic distributions of $M_b$, $R_d$, and $\sigma_0$, and sweeps each of four bias models across literature-plausible parameter ranges. Output reports recovered BTFR zero-point offsets per model.

### §7.5 CMB leakage-coupling bound

```bash
python scripts/cmb_leakage.py \
    --planck-precision 0.005 \
    --output results/cmb_leakage_bound.csv
```

Evaluates the leakage ansatz $g_\text{eff} = g_N + \varepsilon \sqrt{g_N \cdot a_0(z)}$ at the sound-horizon and sub-horizon BAO scales at $z = 1090$, propagates through the Sachs-Wolfe relation, and reports the upper bound on $\varepsilon$ from the Planck first-peak amplitude precision.

### §7.7 Labbé candidate speedup factors

```bash
python scripts/jwst_speedup.py \
    --candidates labbe_2023 \
    --output results/labbe_candidate_speedup.csv
```

Computes $E(z)^{1/4}$ for each of the six central-value massive candidates from the Labbé sample at their photometric redshifts, reporting the multiplicative reduction in required star-formation efficiency under (7.3) of the paper.

---

## Data sources and provenance

| Dataset | Source | Reference |
|---|---|---|
| SPARC rotation curves | [astroweb.cwru.edu/SPARC](http://astroweb.cwru.edu/SPARC/) | Lelli, McGaugh, Schombert, *Astron. J.* **152**, 157 (2016) |
| Übler 2017 KMOS3D BTFR | Tabulated from publication | Übler et al., *Astrophys. J.* **842**, 121 (2017) |
| Labbé 2023 candidate sample | [Public Nature catalog](https://www.nature.com/articles/s41586-023-05786-2) | Labbé et al., *Nature* **616**, 266 (2023) |
| Wisnioski 2015 σ_0 medians | KMOS3D survey paper | Wisnioski et al., *Astrophys. J.* **799**, 209 (2015) |
| Planck 2018 cosmology | [pla.esac.esa.int](https://pla.esac.esa.int/) | Planck Collaboration VI, *Astron. Astrophys.* **641**, A6 (2020) |

The files under `data/` are formatted derivatives of the public sources above, repackaged for direct loading by the analysis scripts. No proprietary data is included.

---

## Reproducibility notes

- **No fitting**: the predictions in this paper are deterministic given the framework's well assignments at $z = 0$ ($\Theta_{a_0} = 13/120$, $\Theta_H = 34/120$) and the standard ΛCDM cosmology (Planck 2018). There are no MCMC chains to converge.
- **Random seeds**: where Monte Carlo is used (the §7.2 Übler forward model), every script accepts a `--seed` argument; default is 42 for deterministic mock samples.
- **Numerical accuracy**: the cosmic-age integral uses SciPy's `quad` with `epsrel=1e-8`; NFW enclosed-mass evaluation is closed-form.
- **Platform**: tested on macOS 14.x and Ubuntu 22.04 with Python 3.11. No GPU required.

---

## License

This repository is released under the MIT License. See `LICENSE` for the full text.

---

## Contact

- Author: B. Shatto, bshatto.pe@gmail.com
- Issues and questions: please open a [GitHub issue](https://github.com/dmobius3/a0z/issues).

For broader work on Mode Identity Theory, see the framework repository at [github.com/dmobius3/mode-identity-theory](https://github.com/dmobius3/mode-identity-theory). For the companion paper on the surface sector (Λcos), see [github.com/dmobius3/lambdacos](https://github.com/dmobius3/lambdacos).
