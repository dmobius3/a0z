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
the bundled location.

## Submission to Foundations of Physics

The current `\documentclass{article}` is for review/draft compilation
only. FoP requires Springer's `svjour3` class at submission time.
`svjour3.cls` is not on CTAN; obtain Springer's submission package
from:

> [https://www.springer.com/journal/10701/submission-guidelines](https://www.springer.com/journal/10701/submission-guidelines)

(look for "LaTeX Macro Package" or "Author Resources"). The bundle
contains `svjour3.cls`, `svglov3.clo`, and `spbasic.bst`. Drop them
into this directory.

### Swap steps

1. In `a0z-paper.tex`, comment out the `\documentclass{article}` line
   and uncomment the `\documentclass[smallextended]{svjour3}` line
   directly below it.

2. Remove `\usepackage{authblk}` from the preamble (svjour3 has its
   own author macros).

3. Replace the current title/author block

    ```latex
    \title{Epoch-Dependent Acceleration Scale from Bounded Topology:\\
           Predictions for High-Redshift Galactic Dynamics}
    \author[1]{B.~Shatto\thanks{...}}
    \affil[1]{Independent Researcher, ...}
    ```

   with the svjour3-flavoured form:

    ```latex
    \title{Epoch-Dependent Acceleration Scale from Bounded Topology:
           Predictions for High-Redshift Galactic Dynamics}
    \author{B. Shatto}
    \institute{B. Shatto \at
        Independent Researcher, St. Petersburg, FL, USA \\
        \email{bshatto.pe@gmail.com} \\
        ORCID: 0009-0007-4419-1311}
    ```

4. Switch the bibliography from inline `thebibliography` back to
   bibtex with Springer's BibTeX style:

    ```latex
    \bibliographystyle{spbasic}
    \bibliography{references}
    ```

   (or keep `thebibliography` inline; svjour3 accepts either).

5. Recompile (`pdflatex`, `bibtex`, `pdflatex`, `pdflatex` if you
   switched to bibtex).

The article-class draft and the svjour3 submission produce the same
content; only the layout/formatting and the title-block macros differ.
