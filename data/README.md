# Observational data

This directory holds formatted derivatives of public observational data products used as inputs to the analysis. No proprietary data is redistributed; refer to the original publications and data archives for license terms.

## Files

| File | Source | Reference |
|---|---|---|
| `lelli_2016_sparc_summary.csv` | [astroweb.cwru.edu/SPARC](http://astroweb.cwru.edu/SPARC/) | Lelli, McGaugh, Schombert, *Astron. J.* **152**, 157 (2016) |
| `ubler_2017_btfr.csv` | Tabulated from publication | Übler et al., *Astrophys. J.* **842**, 121 (2017) |
| `labbe_2023_candidates.csv` | [Nature 616, 266 supplement](https://www.nature.com/articles/s41586-023-05786-2) | Labbé et al., *Nature* **616**, 266 (2023) |

## Schema

### `lelli_2016_sparc_summary.csv`

Local BTFR baseline used as the $z = 0$ reference for §7.2. Columns:

- `galaxy`: SPARC galaxy ID
- `M_b`: baryonic mass [$M_\odot$]
- `v_flat`: asymptotic flat rotation velocity [km/s]
- `v_flat_err`: 1σ uncertainty [km/s]
- `quality`: SPARC quality flag (1, 2, 3)

### `ubler_2017_btfr.csv`

KMOS3D BTFR zero-points at two redshift bins. Columns:

- `z_bin`: nominal redshift (0.9 or 2.3)
- `b`: fitted zero-point ($\log_{10} M_b$ at $\log_{10} v = 2.4$)
- `b_err`: 1σ uncertainty
- `delta_b_obs`: zero-point offset relative to local Lelli baseline
- `n_galaxies`: sample size in this bin

### `labbe_2023_candidates.csv`

The six central-value massive candidates ($\log_{10}(M_\star/M_\odot) > 10$) used for the §7.7 collapse-time speedup analysis. Columns:

- `id`: Labbé catalog ID
- `z_phot`: photometric redshift
- `log_mstar`: $\log_{10}(M_\star/M_\odot)$
- `epsilon_sf_lcdm`: required star-formation efficiency under standard ΛCDM (when reported)

## Usage

Loading scripts in `../scripts/` accept these files via the `--data` argument or default to looking in this directory. Each script's docstring documents the expected schema.
