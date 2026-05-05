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
bibtex   a0z-paper
pdflatex a0z-paper
pdflatex a0z-paper
```

or, equivalently:

```bash
latexmk -pdf a0z-paper.tex
```

Figures are pulled from `../figures/` via `\graphicspath`. Bibliography
is in `references.bib`; bibstyle is `unsrtnat` (citation order matches
the original [1]-[25] numbering of the markdown source).

## Files

| File             | Purpose                                          |
|------------------|--------------------------------------------------|
| `a0z-paper.tex`  | Main LaTeX source                                |
| `references.bib` | BibTeX bibliography (25 entries)                 |
| `README.md`      | This file                                        |

## Notes for arXiv / journal upload

For arXiv, copy `figure{1,2,3,4}.pdf` from `../figures/` into this
directory before tarballing, or adjust `\graphicspath` to point at
the bundled location. The class is plain `article.cls`; Springer
production reformats to their journal style.
