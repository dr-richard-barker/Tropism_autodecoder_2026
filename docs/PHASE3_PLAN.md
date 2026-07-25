# Phase 3 scope — cross-species functional conservation of the spaceflight/tropism signature

Companion to `WEB_TOOL_PLAN.md`. Goal, per the request: for genes altered in the Arabidopsis
tropism/spaceflight signature, reveal which **orthologs are conserved across species** and which
**protein functional regions/functions** are being altered — moving from "which genes overlap" to
"which conserved protein functions are perturbed across kingdoms."

**Decisions (locked):**
- **Functional layer = both** — GO molecular-function/pathway enrichment first, then Pfam/InterPro protein
  **domains** (the literal "functional regions within conserved regions"). See §3 options B + A.
- **Deliverable = precompute + interactive tool view** — static JSON artifacts + a new Phase 3 panel in the
  web tool (§4).

## 1. What already exists (your `OSDR_X-species` bundle — reuse, don't rebuild)
Grounded inventory of `zenodo_repo/`:
- **Orthology hub:** `unified_orthology_matrix.csv` (human Ensembl ↔ species gene, source `babelgene`),
  `arabidopsis_to_human_orthologs.csv` (12,828 AGI → OrthoDB group `og_id` + human Ensembl).
- **Cross-species DEG overlap:** `Table_S2_all_DEGs_with_orthologs.csv` (per-species DEGs on a human-ortholog
  hub) and `Table_S9_cross_species_DEGs.csv` (`human_ortholog, n_species, species_list, directions,
  conservation`) — the **orthology overlap** already computed, incl. direction agreement.
- **Functional layer already done — *location*:** `Table_S3/S4_GOCC_enrichment*`,
  `Table_S5_organelle_enrichment_summary.csv`, `fig4_organelle_enrichment_heatmap`, `fig3_orthology_upset`.
  This tells you **where in the cell** (organelle / GO-Cellular-Component) conserved changes concentrate.
- Species covered: *Arabidopsis thaliana, Homo sapiens, Mus musculus, C. elegans, Drosophila melanogaster*.

**Gap = exactly your ask:** there is **no protein-domain / molecular-function ("what the protein does")
layer**. The existing analysis answers *where*; you want *what functional regions/functions*.

## 2. Key finding that shapes the scope (grounded, not assumed)
- Only **6 / 17** curated tropism markers have a human ortholog; but **~34 % (9,060 / 26,556)** of the
  spaceflight meta-analysis genes do.
- Interpretation: the tropism **receptors/machinery are largely plant-specific** (PINs, phototropins, MIZ1
  have no human ortholog), so the cross-species conservation story is **not** about the tropism sensors —
  it is about the **downstream cellular functions** (stress response, translation, redox, membrane
  transport, proteostasis) that spaceflight perturbs in common across kingdoms. Phase 3 should be framed
  and titled around that, or it will look empty at the receptor level.

## 3. The "functional regions" layer — options (locked: A + B, B first)
| Option | What it shows | Data needed | Effort |
|---|---|---|---|
| **A. Protein domains (Pfam/InterPro)** — *most literal match* | For conserved ortholog groups, the shared **domain families** (kinase, RRM, P-loop NTPase, cytochrome…) and which are over-represented among altered genes | Pfam/InterPro annotations per gene (Ensembl BioMart for human/mouse; InterPro/Araport for Arabidopsis); OrthoDB `og_id` already links the groups | **Medium–High** (new annotation pull + domain enrichment) |
| **B. GO Molecular Function / pathway** | What the conserved altered proteins **do** (catalytic/binding/transporter activity; KEGG/Reactome) | GO-MF + pathway enrichment on the conserved DEG sets — mirrors the existing GOCC pipeline | **Low–Medium** (clusterProfiler, same as the done GOCC) |
| **C. Organelle / GO-CC** | *where* in the cell (already delivered) | — | Done |

**Locked: A + B.** Do **B first** (GO-MF/pathway — fast, reuses the GOCC pipeline) to get the "what
functions are altered" story immediately, then **layer A** (Pfam/InterPro domains) for the sequence-level
"functional regions within conserved regions" detail. Option C (organelle/GO-CC) is already delivered and
is surfaced alongside.

## 4. Proposed architecture (precompute → static → browser)
GitHub Pages is static, so the analysis is **precomputed** and the browser renders/filters it.

**Precompute (Python/R, one-off):**
1. Bridge Arabidopsis signature genes → OrthoDB group (`og_id`) → per-species orthologs (existing tables).
2. Build the conserved-DEG set with per-species direction (from `Table_S9`), tagged conserved / divergent.
3. **Functional annotation** of the conserved set: GO-MF + pathway enrichment (B), and/or Pfam/InterPro
   domain enrichment (A), split by direction (up/down) and species.
4. Emit `docs/assets/ortho/` static JSON: `orthomap.json` (AGI→group→species genes), `conserved_degs.json`,
   `function_enrichment.json` (domain/MF terms × species × direction), plus a species/ortholog index.

**Browser (new Phase 3 view in the tool):**
- **Cross-species input** (original Phase 3): upload non-Arabidopsis data → map to AGI via the ortholog
  table → score with the existing Phase 1/2 engine; report mapping coverage + one-to-many flags.
- **Functional-conservation panel:** for the signature (or the user's altered genes) show (i) an
  orthology-overlap summary (UpSet-style / conservation bars), (ii) a **domain/MF enrichment heatmap**
  (function × species, RdBu, red = up / blue = down), (iii) a drill-down table: conserved gene → shared
  domains → function → per-species direction. Accessible, RdBu, downloadable — same design language as the tool.

## 5. Data sources
- Orthology: your `unified_orthology_matrix.csv` + `arabidopsis_to_human_orthologs.csv` (OrthoDB-based).
- Domains (Option A): Ensembl BioMart `pfam` / `interpro` attributes (human, mouse, fly, worm);
  Araport11/InterPro for Arabidopsis. Ortholog group = OrthoDB `og_id`.
- Function (Option B): GO (org.At.tair.db + org.Hs.eg.db …), KEGG/Reactome via clusterProfiler.

## 6. Staged milestones
- **M1 (data bridge):** signature → ortholog groups → conserved-DEG set; ship `orthomap.json` +
  `conserved_degs.json`; add the cross-species upload/scoring path. *(low risk, mostly existing data)*
- **M2 (function layer — Option B):** GO-MF + pathway enrichment on the conserved set → `function_enrichment.json`
  + the enrichment heatmap panel. *(the "what functions are altered" deliverable)*
- **M3 (domain layer — Option A):** Pfam/InterPro domain annotation + enrichment; add the domain drill-down.
  *(the "functional regions within conserved regions" deliverable)*
- **M4:** polish, legends, downloads; update the roadmap (Phase 3 → live).

## 7. Honest caveats / risks
- Framing: at the tropism-receptor level the overlap is sparse (6/17) — Phase 3 must be presented as
  *conserved cellular functions perturbed by spaceflight*, not conserved tropism sensors.
- Ortholog many-to-many: OrthoDB groups + babelgene give 1-to-many maps; enrichment must handle group-level
  aggregation and report ambiguity, not silently pick one ortholog.
- Domain annotation (Option A) is a real data pull (BioMart/InterPro) and adds weight to `docs/assets/`.
- This reuses a **separate** manuscript's data (`OSDR_X-species`) — confirm licensing/attribution before
  bundling its tables into this repo, or link to its Zenodo record instead of copying.
