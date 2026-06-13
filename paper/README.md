# a0z-paper.tex

LaTeX source for the preprint

> **B. Shatto, "Epoch-Dependent Acceleration Scale from Bounded Topology:
> Predictions for High-Redshift Galactic Dynamics" (2026).**

Source of truth for the prose remains the markdown copy at
`mode-identity-theory/files/working/files/a0-evolution-paper.md`; this
directory holds the typeset LaTeX rendering for open-archive deposit.

## Build

```bash
cd paper
pdflatex a0z-paper
pdflatex a0z-paper
```

Two passes are sufficient. The bibliography is inlined as a
`thebibliography` environment, so no bibtex pass is required. Figures
are pulled from `../figures/` via `\graphicspath`; `references.bib` is
kept alongside for optional bibtex regeneration.

## Class

The preprint builds with the standard `article` class
(`\documentclass[11pt,a4paper]{article}`) for a clean, journal-neutral
typeset (18 pages).

The Springer `svjour3` bundle (`svjour3.cls`, `svglov3.clo`,
`sp{basic,phys,mpsci}.bst`) is retained in this directory from an
earlier formatting pass, in case of a future Springer-journal
resubmission. It is **not** used by the current build.

## Files

| File | Purpose |
|------|---------|
| `a0z-paper.tex` | Main LaTeX source (article class) |
| `references.bib` | BibTeX bibliography (25 entries, optional) |
| `svjour3.cls`, `svglov3.clo`, `sp*.bst` | Springer bundle, retained but unused by current build |
| `README.md` | This file |
| `.gitignore` | LaTeX intermediates |

## Building a source bundle for deposit

Most archives (SSRN, Zenodo) accept the compiled `a0z-paper.pdf`
directly. For an archive that compiles from source, copy
`figure{1,2,3,4}.pdf` from `../figures/` into this directory (or adjust
`\graphicspath`) before tarballing the `.tex` + `references.bib` +
figures.
