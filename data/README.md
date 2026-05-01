# Observational data

This directory holds formatted derivatives of public observational data products used as inputs to the analysis. No proprietary data is redistributed; refer to the original publications and data archives for license terms.

## Files

| File | Source | Reference |
|---|---|---|
| `labbe_2023_candidates.csv` | [Nature 616, 266 supplement](https://www.nature.com/articles/s41586-023-05786-2) | Labbé et al., *Nature* **616**, 266 (2023) |

## Schema

### `labbe_2023_candidates.csv`

The six central-value massive candidates ($\log_{10}(M_\star/M_\odot) > 10$) used for the §3.4 collapse-time speedup analysis. Columns:

- `id`: Labbé catalog ID
- `z_phot`: photometric redshift
- `log_mstar`: $\log_{10}(M_\star/M_\odot)$

The §3.4 prediction $\epsilon_\text{SF}^\text{MIT} = \epsilon_\text{SF}^\text{ΛCDM} / E(z)^{1/4}$ is computed by `jwst_speedup.py` directly from `z_phot`; it does not consume a per-candidate $\epsilon_\text{SF}^\text{ΛCDM}$ column.

## Other observational inputs

The SPARC (Lelli et al. 2016), Übler+ 2017 KMOS3D BTFR, Wisnioski 2015 σ_0 medians, and Planck 2018 cosmology values used elsewhere in the pipeline are taken directly from the cited publications and hardcoded in the relevant scripts (`framework.py`, `ubler_sigma_tension.py`, `cosmology.py`). They are not redistributed here as separate CSV files.
