# LaTeX manuscript (npj Microgravity / Springer Nature style)

Assembles the manuscript in the official **Springer Nature LaTeX template**
(`sn-jnl` class) — the format npj Microgravity accepts and typesets from.

```
latex/
├── main.tex          # assembled manuscript (sn-jnl class, sn-nature refs)
├── references.bib    # references transcribed from the manuscript
├── figures/          # the 13 figures (PNG)
└── README.md         # this file
```

## How to compile

`sn-jnl.cls` / `sn-nature.bst` ship with Springer Nature's official template
(not vendored here).

- **Overleaf (recommended):** new project from the **"Springer Nature Article
  Template (sn-jnl)"** → replace its `main.tex` with this one, upload
  `references.bib` and `figures/`, compile with **pdfLaTeX**.
- **Local:** place `sn-jnl.cls` + `sn-nature.bst` here, then
  `pdflatex main` → `bibtex main` → `pdflatex main` → `pdflatex main`.

## Status / TODO before submission

- [x] **Compile-tested** — builds with latexmk (pdflatex + BibTeX) without errors or BibTeX warnings (28 Oct 2026); CI also builds it.
- [x] **Authorship** — confirmed and corrected (removed AI tool "Phylo Biomni" as author).
- [x] **References** — checked against CrossRef and completed with volume/pages/DOI (28 Oct 2026); four wrong entries corrected.
- [ ] **Figures** — repo PNGs. For final submission npj prefers vector (PDF/EPS)
      or ≥300 dpi; swap files in `figures/` and paths resolve unchanged.

## Source

Ported from `../manuscript.md`. Body text, figures, and references are the
author's own content — nothing was invented.
