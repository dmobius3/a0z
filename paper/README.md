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
regenerating the bibliography via bibtex if desired (replace the
inline `\begin{thebibliography}...\end{thebibliography}` with
`\bibliographystyle{unsrtnat}` plus `\bibliography{references}` and
add a bibtex pass to the build).

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
