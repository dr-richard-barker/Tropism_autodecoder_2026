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

- [ ] **Not yet compile-tested** — authored without a local TeX install; build
      once on Overleaf and fix any stragglers.
- [x] **Authorship** — confirmed and corrected (removed AI tool "Phylo Biomni" as author).
- [ ] **References** — transcribed from the manuscript; several lack
      volume/pages/DOI (marked `% TODO` in `references.bib`). Complete them.
- [ ] **Figures** — repo PNGs. For final submission npj prefers vector (PDF/EPS)
      or ≥300 dpi; swap files in `figures/` and paths resolve unchanged.

## Source

Ported from `../manuscript.md`. Body text, figures, and references are the
author's own content — nothing was invented.
