# Arabidopsis Spaceflight Tropism Recognition Pipeline — Implementation Plan

## Summary

Build a reproducible, FAIR-compliant pipeline that (1) systematically ingests and harmonizes all Arabidopsis thaliana spaceflight transcriptomics from NASA OSDR + complementary GEO, (2) trains a custom variational **stimulus auto-decoder** on the Salk Arabidopsis Developmental Atlas (GSE226097, 432,919 nuclei) to learn cell-type- and stimulus-specific latent codes, (3) deconvolves every bulk spaceflight sample onto the atlas to recover cell-type-specific tropism response scores, (4) trains a tropism-condition classifier spanning gravitropism, phototropism, thigmotropism/mechanotropism, hydrotropism, chemotropism, and oxytropism, and (5) renders three complementary visualization tracks: KEGG pathway figures (ggkegg), machine-interpretable SBGN-ML maps (SBGNview), and **tissue-anatomy heatmaps (ggPlantMap)** that project the cell-type-specific stress-decoding scores onto Arabidopsis organ/cell-type maps to show *which cell types are perceiving altered tropic signalling*. Output is a Zenodo-ready repository with Manuscript/Code folders, DataCite metadata.json, and a machine-readable README.

### Honest scope note (data asymmetry)
Tropism transcriptomics is uneven. Gravitropism and phototropism are well-covered by spaceflight + ground datasets. Hydrotropism has dedicated transcriptomics (brassinosteroid/ecotype study, PubMed 29439211). Thigmotropism/mechanotropism is represented indirectly via mechanosensitive-channel and FER/PIF3 root-penetration datasets (GSE225299) rather than a canonical "thigmotropism" series. Chemotropism and oxytropism have **no dedicated Arabidopsis transcriptomics** in GEO/OSDR. The pipeline will: (a) include all available data per tropism, (b) flag tropisms with insufficient data as "reference-signature-only" (curated marker gene sets from literature, scored via GSVA/aucell rather than trained from data), and (c) report coverage transparently in the manuscript. No data will be fabricated.

---

## Key decisions (locked from clarification)

| Axis | Decision |
|---|---|
| Recognition target | Both — cell-type-specific deconvolution scores AND a sample-level tropism-condition classifier |
| Data scope | Systematic pull of all Arabidopsis spaceflight studies from OSDR + complementary GEO tropism series |
| Auto-decoder depth | Custom variational auto-decoder trained on the Salk atlas (encoder→stimulus latent codes→decoder) |
| Deliverable scope | Full end-to-end run: download everything, train on full 400k+ atlas, deconvolve + visualize every study, generate all figures |
| Tropism coverage | Gravitropism, phototropism, thigmotropism/mechanotropism, hydrotropism, chemotropism, oxytropism — with honest coverage flags |

---

## Verified external resources

- **Salk Arabidopsis Developmental Atlas** = GEO **GSE226097** (Lee, Nobori, Illouz-Eliaz et al., *Nature Plants* 2025). 432,919 nuclei, 10 developmental stages, 29 supplementary Seurat `.rds` files including `GSE226097_global_integration_221009.rds` (~3.2 GB) and per-organ objects. MERFISH spatial datasets also available. [12, 16, 50, 51]
- **NASA OSDR Biological Data API** — REST + query interface, ISA-Tab model. Confirmed live: `https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/` returns 62 Arabidopsis spaceflight studies, 2,204 samples (874 flight, 400 ground control), 20 RNA-seq + 6 microarray studies with assay metadata. [21, 23, 25]
- **PhysioSpace / Plant PhysioSpace** — Lenz et al. 2013 (PLOS ONE) and Esfahani et al. 2021 (*Plant Physiology*); spherical transform + Wilcoxon PhysioScores, cross-platform/cross-species. Used as a baseline comparator, not the primary method. [14, 15]
- **ggkegg** (Bioconductor, Sato) — `pathway()` fetches KGML → `tbl_graph`, integrates with clusterProfiler; produces KEGG-native maps via ggraph. Does NOT emit SBGN. [1, 3, 8]
- **SBGNview** (Bioconductor) + **KEGGtranslator** — render SBGN-ML to SVG, support KEGG/Reactome/MetaCyc; KEGGtranslator converts KGML→SBGN-ML. This is the SBGN-compliant stack. [41, 43, 45, 48]
- **ggPlantmap** (Jo & Kajala 2024, *J Exp Bot*, doi:10.1093/jxb/erae043; GitHub `leonardojo/ggPlantmap`) — open-source R package for eFP-like quantitative heatmaps projected onto plant tissue/cell-type maps via ggplot2. Pre-loaded Arabidopsis maps directly relevant to tropism tissues: `ggPm.At.roottip.crosssection`, `ggPm.At.roottip.longitudinal`, `ggPm.At.rootelong.longitudinal`, `ggPm.At.rootmatur.crosssection`, `ggPm.At.3weekrosette.topview`, `ggPm.At.leaf.crosssection`, `ggPm.At.leafepidermis.topview`, `ggPm.At.shootapex.longitudinal`, `ggPm.At.inflorescencestem.crosssection`, `ggPm.At.seed.devseries`, `ggPm.At.earlyembryogenesis.devseries`. Workflow: `ggPlantmap.merge(map, quant, id.x="ROI.name")` joins per-cell-type scores to tissue ROIs, `ggPlantmap.plot()` renders the heatmap onto the anatomy image; custom maps can be built from plant images via XML for organs not pre-loaded (e.g. apical hook, hypocotyl — key gravitropism/phototropism organs). [54, 55, 57, 61]
- **CIBERSORTx** — validated for Arabidopsis bulk deconvolution (Vong et al. 2024, *Plant Physiology*; Moreno 2024). Used as a baseline comparator against the custom auto-decoder. [13, 18]

## Environment state (verified)

- Default machine `worker-0`: 1 CPU, 16 GB RAM, **no GPU**. R 4.x + Python 3.11 in `/opt/conda` and `/workspace/.venv`.
- Installed R: DESeq2, limma, clusterProfiler, GEOquery, SingleCellExperiment, Seurat, ComplexHeatmap, KEGGREST, fgsea, ggplot2, svglite, tidyverse, data.table, BiocManager.
- **Missing R (will install)**: ggkegg, SBGNview, tidygraph, ggraph, org.At.tair.db, pathview, SingleR, scran, scater, reactome.db, **ggPlantmap** (from GitHub `leonardojo/ggPlantmap` since it is not on CRAN/Bioc).
- Installed Python: scanpy 1.11.4, anndata 0.12.1, torch 2.7.1+cu126 (CPU only here), pytorch_lightning 2.5.2.
- Nextflow available at `/usr/local/bin/nextflow`. Docker present.
- No GPU on default machine → auto-decoder training needs a GPU machine via `ManageMachine`.

---

## Compute & resource estimate

| Stage | Input | RAM | Disk | Runtime | Target |
|---|---|---|---|---|---|
| OSDR + GEO metadata + raw download | 26 transcriptomics studies (~1,300 samples) | 4 GB | ~150–300 GB raw FASTQ | 2–6 h (network-bound) | right-sized sandbox (8 CPU/32 GB) |
| Salk atlas download + load | ~3.2 GB global Seurat `.rds` (432k nuclei) | 24–32 GB to load | 8 GB | 20–40 min download + load | 32 GB machine |
| Atlas preprocessing (normalize, HVG, cell-type labels) | 432k × ~25k genes | 16–24 GB | 4 GB | 30–60 min | 32 GB machine |
| Auto-decoder training (VAE, latent stimulus codes) | 432k nuclei subset to HVG ~3–5k genes | 8–16 GB GPU RAM | 2 GB | 2–6 h on 1 GPU | **GPU machine** (ManageMachine) |
| Per-study bulk DE + deconvolution | 26 studies × ~50 samples × ~25k genes | 8 GB | 1 GB | 5–15 min/study → ~3–6 h total | fan-out across 2 machines |
| Tropism classifier training | deconvolution scores × tropism labels | 4 GB | <1 GB | 20–40 min | default machine |
| ggkegg + SBGNview + ggPlantMap figures | DE tables + KEGG/Reactome pathways + cell-type tropism scores | 4 GB | 2 GB figures | 1.5–2.5 h | default machine |
| Repo assembly + metadata.json | all outputs | 2 GB | <1 GB | 15 min | default machine |

**Execution target decision**: Use `ManageMachine` to provision (a) one 8-CPU/32-GB machine for atlas download/load/preprocessing and OSDR bulk download, and (b) one GPU machine for auto-decoder training. Fan out per-study deconvolution across the two allowed machines. Total wall-clock estimate: **8–14 h** end-to-end with chunked checkpoints to `/mnt/shared-workspace/`. This fits within the 24 h sandbox cap if chunked; auto-decoder training is the single longest step and will checkpoint model weights to shared-workspace.

---

## Pipeline architecture (behavior-level)

### Subsystem 1 — Data acquisition & harmonization (`Code/01_data_acquisition/`)
- **OSDR puller** (`pull_osdr.py`): query `v2/query/metadata/` for `study.characteristics.organism=Arabidopsis thaliana AND study.factor value.spaceflight`; enumerate the 26 transcriptomics studies (20 RNA-seq + 6 microarray); pull per-study ISA-Tab metadata + raw/processed expression files via the OSDR file API; cache to `/mnt/shared-workspace/raw_osdr/`.
- **GEO puller** (`pull_geo.py`): use `GEOquery` to fetch complementary tropism series (gravitropism GSE3847/GSE8300 family, phototropism GSE143760, hydrotropism ecotype study, mechanotropism GSE225299); cache to `/mnt/shared-workspace/raw_geo/`.
- **Atlas puller** (`pull_atlas.py`): download `GSE226097_global_integration_221009.rds` + per-organ Seurat objects from GEO FTP; cache to `/mnt/shared-workspace/salk_atlas/`.
- **Metadata schema** (`metadata_schema.json` + `harmonize_metadata.py`): a unified sample-level schema linking OSDR sample IDs → spaceflight telemetry (factor value.spaceflight, hardware, mission, growth temperature, light regime, gravity condition, tissue, ecotype, developmental stage) and GEO sample IDs → tropism condition. Output: `harmonized_metadata.tsv` with columns `sample_id, source_repo, source_accession, organism, tissue, ecotype, developmental_stage, spaceflight_condition, gravity, light_regime, tropism_type, hardware, mission, assay_type, raw_file, processed_file`.

### Subsystem 2 — Atlas preprocessing & reference building (`Code/02_atlas_reference/`)
- **Load & convert** (`prepare_atlas.py`): load the global Seurat `.rds` via `Seurat`/`scholar` or convert to AnnData with `scanpy`; retain 432,919 nuclei; map cluster labels (183 major + 653 subclusters) to a consolidated cell-type ontology.
- **Normalize & HVG** (`scanpy`): SCTransform or log-normalize; select 3,000–5,000 HVGs; store as `atlas_reference.h5ad` in `/workspace/` then copy to `/mnt/shared-workspace/`.
- **Cell-type signature matrix** (`build_signatures.py`): per-cell-type mean expression + marker genes (top 50/cell type); export `cell_type_signatures.csv` for CIBERSORTx baseline and PhysioSpace scoring.

### Subsystem 3 — Custom stimulus auto-decoder (`Code/03_auto_decoder/`)
- **Architecture** (`model.py`, PyTorch + pytorch_lightning): a conditional VAE / auto-decoder where the encoder maps a cell's expression vector → latent `z` (stimulus code, dim 16–32), and the decoder reconstructs expression conditioned on `(z, cell_type_code, developmental_stage_code)`. Trained on the atlas with cell-type and stage as conditioning labels. Stimulus codes are regularized toward Gaussian prior (KL) so that held-out stimuli project to interpretable regions.
- **Training** (`train.py`): train on atlas HVG matrix; checkpoint weights + latent embeddings every epoch to `/mnt/shared-workspace/autodecoder/`; log reconstruction loss + KL + cell-type classification accuracy (auxiliary head). Target: 50–100 epochs, 2–6 h on GPU.
- **Bulk projection** (`project_bulk.py`): for each bulk sample, optimize a per-sample stimulus latent code (and cell-type mixture weights) that minimizes reconstruction error between the bulk profile and the decoder output weighted by estimated cell-type proportions; output `stimulus_latent_per_sample.tsv` + `cell_type_proportions.tsv` + `cell_type_specific_tropism_scores.tsv` (PhysioScore-style signed log10 p per cell type per sample).
- **Baseline comparators** (`baselines/`): CIBERSORTx (signature matrix from atlas) and Plant PhysioSpace (cross-platform stimulus scores) run on the same bulk samples; report concordance (Spearman, RMSE) vs. the auto-decoder.

### Subsystem 4 — Differential expression & tropism classifier (`Code/04_de_and_classifier/`)
- **Per-study DE** (`de_per_study.py`): DESeq2 for RNA-seq studies, limma for microarray; contrast = spaceflight vs. ground control within each OSDR study; output per-study DE tables (`de_<OSD-id>.tssv`) with gene, log2FC, padj.
- **Meta-analysis** (`meta_de.py`): combine across studies with `metafor` (random-effects) on log2FC; output `meta_de_gravitropism.tsv` etc. partitioned by tropism type using metadata.
- **Tropism classifier** (`tropism_classifier.py`): features = auto-decoder stimulus latent codes + cell-type proportions + per-cell-type tropism scores; labels = tropism condition from harmonized metadata; model = elastic-net logistic regression with **nested 5-fold CV** (feature selection inside folds), reporting mean outer-fold macro-F1 + per-class F1 + confusion matrix. For tropisms with insufficient labeled samples (chemotropism, oxytropism), use curated literature marker signatures scored via GSVA/AUCell as a semi-supervised fallback, flagged as "reference-signature-only."

### Subsystem 5 — Systems biology & tissue-anatomy visualization (`Code/05_visualization/`)

Three complementary tracks, each producing `.svg` (primary, editable) + `.png` (300 dpi) to `/mnt/results/figures/`:

**Track A — KEGG-native pathway maps (ggkegg)** (`ggkegg_figures.R`): for each tropism condition, fetch relevant KEGG pathways (ath04075 = plant hormone signal transduction, ath00040 = pentose phosphate, ath00900 = terpenoid biosynthesis, ath02010 = ABC transporters, ath00071 = fatty acid degradation, etc.); overlay log2FC from meta-DE onto pathway nodes via `pathway()` + `tidygraph` + `ggraph`; export with `ggkeggsave`. Save to `/mnt/results/figures/ggkegg/`.

**Track B — SBGN-ML machine-interpretable maps (SBGNview)** (`sbgn_figures.R`): use `SBGNview` to render SBGN-ML maps for the same pathways (KEGG→SBGN via KEGGtranslator or SBGNview's KEGG collection); overlay expression + predicted metabolites/cofactors (from KEGGREST compound mapping) as glyph colors; export `.svg` (machine-interpretable, validates against SBGN schema) + `.png`. Save to `/mnt/results/figures/sbgn/`.

**Track C — Tissue-anatomy tropism-decoding heatmaps (ggPlantMap)** (`ggplantmap_figures.R`): project the **cell-type-specific tropism/stress-decoding scores** from Subsystem 3 onto Arabidopsis tissue maps to show *which cell types are perceiving altered tropic signalling*. This complements the pathway maps (Tracks A/B show molecular circuits; Track C shows where in the plant the decoding lands anatomically).
- **Cell-type → ROI mapping** (`build_celltype_to_roi.py`): build a controlled mapping from the Salk atlas cell-type labels (183 major clusters → consolidated cell types) to ggPlantMap `ROI.name` values per tissue map (e.g. atlas "Columella" → `ggPm.At.roottip.crosssection` ROI "Columella"; atlas "Trichoblast/Atrichoblast" → root epidermis ROIs; atlas "Xylem/Phloem/Pericycle" → vascular ROIs; atlas "Mesophyll/Palisade" → `ggPm.At.leaf.crosssection` ROIs). Unmappable cell types are flagged and either grouped to the nearest anatomical ROI or rendered on a companion UMAP (Subsystem 5 summary figures) rather than forced onto a tissue.
- **Per-condition rendering**: for each (tropism condition × tissue) pair, aggregate the auto-decoder's per-cell-type tropism scores (mean signed PhysioScore across samples in that condition) and render via `ggPlantmap.merge(map, quant, id.x="ROI.name")` + `ggPlantmap.plot()` with a diverging colorblind-friendly palette (blue = downregulated decoding, red = upregulated). Tissues mapped: root tip (cross-section + longitudinal), root elongation/mature zones, rosette top view, leaf cross-section + epidermis, shoot apex, inflorescence stem, seed/embryo series.
- **Custom maps for tropism-key organs** (`build_custom_maps.py` + XML): if the pre-loaded maps lack apical hook and hypocotyl (central to gravitropism/phototropism), build custom ggPlantMap XML maps from published Arabidopsis hypocotyl/hook anatomy images (the Salk atlas itself includes apical-hook spatial profiling, so we can trace ROIs from their published figures) and register them as `ggPm.At.hypocotyl.longitudinal` and `ggPm.At.apicalhook.longitudinal`.
- **Multi-panel composite figures** (`composite_figures.R`): for each tropism, assemble a publication figure combining (top) the ggPlantMap tissue heatmap with cell-type tropism scores, (middle) the ggkegg pathway overlay for the dominant signalling pathway, (bottom) the SBGN-ML map for machine-interpretable export — linked by shared gene/metabolite IDs. Save to `/mnt/results/figures/composite/`.

**Summary figures** (`summary_figures.py`): UMAP of atlas with stimulus latent overlay; heatmap of cell-type-specific tropism scores across conditions (the "flat" companion to the ggPlantMap tissue projection); confusion matrix for the classifier; volcano of meta-DE per tropism; metabolite/cofactor network. All `.svg` + `.png`.

### Subsystem 6 — FAIR repository assembly (`Code/06_repository/`)
- **Directory structure**:
  ```
  arabidopsis_spaceflight_tropism/
  ├── README.md                      # machine-readable manifest (YAML front matter + human text)
  ├── metadata.json                  # DataCite schema 4.x
  ├── environment.yml                # conda env for Python+R
  ├── Dockerfile                     # containerized runtime
  ├── nextflow.config                # Nextflow config
  ├── main.nf                        # Nextflow entrypoint orchestrating 01–06
  ├── Code/
  │   ├── 01_data_acquisition/
  │   ├── 02_atlas_reference/
  │   ├── 03_auto_decoder/
  │   ├── 04_de_and_classifier/
  │   ├── 05_visualization/
  │   │   ├── ggkegg_figures.R
  │   │   ├── sbgn_figures.R
  │   │   ├── ggplantmap_figures.R
  │   │   ├── build_celltype_to_roi.py
  │   │   ├── build_custom_maps.py        # custom ggPlantMap XML for hypocotyl/apical hook
  │   │   ├── composite_figures.R
  │   │   ├── summary_figures.py
  │   │   └── custom_maps/                # custom ggPlantMap XML + source images
  │   ├── 06_repository/
  │   └── tests/                     # pytest + testthat smoke tests
  ├── Manuscript/
  │   ├── 00_main.md                 # npj Microgravity template placeholders
  │   ├── 01_methods.md
  │   ├── 02_results.md
  │   ├── 03_discussion.md
  │   ├── figures/                   # high-res .svg + .png + legends
  │   │   ├── ggkegg/                 # KEGG-native pathway maps
  │   │   ├── sbgn/                   # SBGN-ML machine-interpretable maps
  │   │   ├── ggplantmap/             # tissue-anatomy tropism-decoding heatmaps
  │   │   └── composite/              # multi-panel tropism summary figures
  │   ├── tables/                    # CSV/TSV
  │   └── references.bib
  ├── Data/
  │   ├── harmonized_metadata.tsv
  │   ├── de_results/
  │   ├── deconvolution_results/
  │   ├── tropism_scores/
  │   └── celltype_to_roi_mapping.tsv  # atlas cell types → ggPlantMap ROIs
  └── docs/
      ├── API.md                     # OSDR/GEO/KEGG API usage documented
      ├── SBGN_compliance.md
      └── FAIR_checklist.md
  ```
- **metadata.json** (DataCite 4.x): creators, titles, publisher (Zenodo), publicationYear, resourceType (ComputationalNotebook/Software), subjects (space biology, Arabidopsis, tropism, scRNA-seq deconvolution, SBGN), rights (CC-BY-4.0), relatedIdentifiers (GSE226097, OSDR study accessions, GitHub repo), version, descriptions, formats, sizes, fundingReference (NASA if applicable).
- **README.md**: YAML front matter (`id`, `title`, `version`, `pipeline`, `inputs`, `outputs`, `run_command`, `license`, `doi_placeholder`) followed by human-readable overview, quickstart, methods summary, figure index, and FAIR compliance statement.
- **Reproducibility**: `environment.yml` pins all Python + R packages; `Dockerfile` builds on `bioconductor/bioconductor_docker:RELEASE_3_20` + pytorch; `nextflow.config` defines profiles (`local`, `docker`, `hpc`); `main.nf` chains the six subsystems with checkpointing.

---

## Testing & acceptance criteria

- **Smoke tests** (`Code/tests/`): (1) OSDR query returns ≥60 Arabidopsis spaceflight studies; (2) atlas `.rds` loads and has >400k nuclei; (3) auto-decoder trains and reconstruction loss decreases monotonically over epochs; (4) deconvolution of a held-out atlas pseudobulk recovers known cell-type proportions within 10% MAE; (5) ggkegg renders ≥5 pathway figures without error; (6) SBGNview renders ≥5 SBGN-ML files that validate against the SBGN schema; (7) **ggPlantMap renders ≥3 tissue heatmaps (root tip, leaf, hypocotyl) with cell-type tropism scores mapped to ROIs, and the cell-type→ROI mapping table covers ≥80% of atlas major cell types**; (8) `metadata.json` validates against DataCite 4.x JSON schema; (9) `nextflow run main.nf -profile test` completes on a 2-study subset.
- **Acceptance**: all 26 transcriptomics studies processed; auto-decoder trained with documented loss curves; deconvolution + classifier outputs for every sample; ≥10 ggkegg + ≥10 SBGN + **≥6 ggPlantMap tissue heatmaps (covering root, leaf, shoot apex, hypocotyl/apical hook, stem, seed)** + ≥4 composite multi-panel figures; complete repo with valid metadata.json; README renders both YAML and human sections.

## Assumptions

- OSDR file-download endpoints behave as documented; if a study's raw FASTQ is unavailable, fall back to OSDR-processed count matrices (GeneLab processed-data pipeline GL-DPPD).
- The global Seurat `.rds` loads in R on a 32 GB machine; if memory-bound, fall back to per-organ `.rds` files and integrate with `Seurat::merge()`.
- GPU machine is provisionable via `ManageMachine` within the Free tier (2 machines/session). If GPU unavailable, train auto-decoder on CPU at reduced scale (50k-nucleus subsample, fewer epochs) and flag the reduction.
- Tropism labels for OSDR samples are derived from metadata `study.factor value.spaceflight` + tissue + light regime + hardware; ambiguous samples are labeled "unspecified_tropism" and excluded from classifier training but retained for deconvolution.
- Chemotropism/oxytropism classifier classes are reference-signature-only due to missing transcriptomics; this is reported as a limitation, not hidden.
- **ggPlantMap cell-type → ROI mapping**: the Salk atlas cell-type labels (183 major clusters) will be consolidated to a vocabulary that aligns with ggPlantMap's pre-loaded `ROI.name` sets; clusters without a clear anatomical home (e.g. uncharacterized states) are rendered on the companion UMAP rather than forced onto a tissue. Custom maps for hypocotyl/apical hook will be built from the Salk atlas's own published apical-hook spatial figures (the atlas paper profiles this organ explicitly), so ROI boundaries are traceable to a primary source.

## Next steps after approval

1. Provision machines via `ManageMachine` (32 GB sandbox + GPU sandbox).
2. Install missing R packages (ggkegg, SBGNview, tidygraph, ggraph, org.At.tair.db, pathview, SingleR, scran, scater) and Python deps (scvi-tools optional, pytorch-lightning already present).
3. Execute subsystems 01 → 06 in order, checkpointing to `/mnt/shared-workspace/` between stages and surfacing figures/tables to `/mnt/results/` as they land.
4. Assemble final repo in `/mnt/results/arabidopsis_spaceflight_tropism/` and zip for Zenodo upload.
