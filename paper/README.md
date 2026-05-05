# a0z-paper.tex

LaTeX source for the paper

> **B. Shatto, "Epoch-Dependent Acceleration Scale from Bounded Topology:
> Predictions for High-Redshift Galactic Dynamics" (2026).**

Target journal: Foundations of Physics. Source of truth for the prose
remains the markdown copy at
`mode-identity-theory/files/working/files/a0-evolution-paper.md`; this
directory holds the LaTeX rendering for journal submission.

## Build

```bash
cd paper
pdflatex a0z-paper
pdflatex a0z-paper
```

Two passes are sufficient. The bibliography is inlined as a
`thebibliography` environment so no bibtex pass is required.

Figures are pulled from `../figures/` via `\graphicspath`. The
companion `references.bib` is kept alongside as a convenience for
regenerating the bibliography via bibtex if desired.

## Class

The paper builds with Springer's **`svjour3`** class
(`\documentclass[smallextended]{svjour3}`), the official LaTeX class
for Foundations of Physics submissions. The class file plus its
companions (`svglov3.clo`, `spbasic.bst`, `spphys.bst`, `spmpsci.bst`)
are bundled in this directory, sourced from
[github.com/latextemplates/svjour](https://github.com/latextemplates/svjour).
Upstream copyright is Springer; the files are only valid for
Springer-journal submissions.

There is a commented-out `\documentclass{article}` line directly below
the active svjour3 line for review/draft compilation if you want to
bypass the Springer layout.

## Files

| File             | Purpose                                           |
|------------------|---------------------------------------------------|
| `a0z-paper.tex`  | Main LaTeX source                                 |
| `references.bib` | BibTeX bibliography (25 entries, optional)        |
| `svjour3.cls`    | Springer journal class                            |
| `svglov3.clo`    | Class options file                                |
| `sp{basic,phys,mpsci}.bst` | Springer BibTeX styles                  |
| `README.md`      | This file                                         |
| `.gitignore`     | LaTeX intermediates                               |

## Notes for arXiv / journal upload

For arXiv, copy `figure{1,2,3,4}.pdf` from `../figures/` into this
directory before tarballing, or adjust `\graphicspath` to point at
the bundled location. Include `svjour3.cls`, `svglov3.clo`, and
`sp*.bst` in the submission tarball (arXiv won't have them otherwise).
