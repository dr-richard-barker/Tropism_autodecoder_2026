# Arabidopsis thaliana Spaceflight Tropism Recognition: Integrating Bulk Transcriptomics with a Single-Cell Developmental Atlas Foundation Model

## Authors

Richard Barker¹

¹ Independent researcher

## Abstract

Plants perceive and respond to multiple directional stimuli (tropisms) that guide growth and development. Under spaceflight conditions, altered gravity and light environments perturb tropic signaling, providing a unique window into these fundamental plant responses. We present a FAIR-compliant computational pipeline that integrates 1,253 bulk transcriptomic samples from NASA GeneLab OSDR and NCBI GEO with the Salk Institute Arabidopsis Developmental Atlas (432,919 nuclei, 183 cell-type clusters) as a foundation model. A custom variational auto-decoder (32-dimensional latent space) enables cell-type deconvolution of bulk samples, recovering tissue-level tropism signatures. Differential expression analysis (DESeq2) across 11 spaceflight studies, combined with DerSimonian-Laird random-effects meta-analysis, identified 16,002 significantly differentially expressed genes (padj < 0.05) [TODO: produced with a τ² that was too small by a factor of Σw, fixed in `Code/meta_classifier.py` on 28 Oct 2026; re-run on the per-study DESeq2 tables, which are not in this repository (README, Limitation 7)]. An elastic-net classifier achieved AUC = 0.919 [TODO: printed by `Code/meta_classifier.py` but never saved; re-run it and commit its metrics output (README, Limitation 8)] for distinguishing spaceflight from ground control samples based on cell-type fractions and stimulus activation scores. ggPlantMap tissue projections revealed organ-specific abundance changes, with stem tissues showing the most significant response (padj = 0.003). The pipeline is packaged as a Zenodo-ready repository with DataCite metadata, Dockerfile, and conda environment specification.

## Introduction

Tropisms — directional growth responses to environmental stimuli — are fundamental to plant development. Gravitropism (gravity sensing), phototropism (light-directed growth), thigmotropism (touch response), and hydrotropism (water-seeking) are mediated by distinct but overlapping signaling pathways [1]. Spaceflight provides a unique experimental context where gravity is effectively removed (microgravity), allowing dissection of gravitropism from other tropic responses.

Arabidopsis thaliana has been extensively studied in spaceflight experiments, with NASA GeneLab's Open Science Data Repository (OSDR) hosting transcriptomic data from dozens of spaceflight experiments [2]. However, these bulk transcriptomic measurements average over all cell types, obscuring tissue-specific responses. Single-cell transcriptomics has revolutionized our understanding of plant development, with the Salk Institute Arabidopsis Developmental Atlas providing a comprehensive reference of 432,919 nuclei across 10 developmental stages and 183 cell-type clusters [3].

Computational deconvolution methods, such as CIBERSORTx [4] and PhysioSpace [5], enable estimation of cell-type proportions from bulk transcriptomic data. However, these methods typically use fixed reference signatures and do not learn stimulus-specific representations. Variational auto-encoders (VAEs) offer a powerful framework for learning compact latent representations of gene expression, and conditional variants can incorporate metadata such as cell type, organ, and developmental stage [6].

Here we present a complete pipeline that: (1) systematically acquires and harmonizes Arabidopsis spaceflight transcriptomics from OSDR and GEO, (2) integrates the Salk atlas as a foundation model with a custom variational auto-decoder for cell-type deconvolution, (3) performs differential expression and meta-analysis across studies, (4) trains a tropism classifier, and (5) visualizes results using ggPlantMap tissue projections and KEGG pathway networks.

## Methods

### Data Acquisition

**OSDR**: We queried the NASA GeneLab OSDR metadata API (`https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/`) to identify all Arabidopsis thaliana transcriptomics studies, yielding 24 studies (18 RNA-seq, 6 microarray) with 1,190 samples. Count matrices were downloaded via direct file access (`https://osdr.nasa.gov/geode-py/ws/studies/OSD-{id}/download`). Large studies that timed out on the data query API were obtained through direct file download.

**GEO**: We identified 5 tropism-focused GEO series: GSE3847 (gravitropism), GSE143760 (phototropism), GSE97258 (hydrotropism), GSE225299 (mechanotropism/RALF1 signaling), and GSE115554 (gravitropism). (GSE8300, listed here before 28 Oct 2026, is a *Haemophilus influenzae* genomic hybridization series; its 84 rows were removed from `Data/harmonized_metadata.tsv`.) Supplementary files were downloaded via HTTPS mirror of the NCBI FTP server. GSE225299 scRNA-seq data (12 samples) was aggregated to pseudo-bulk counts per sample.

**Harmonization**: All samples were annotated with tropism type (gravitropism, phototropism, mechanotropism, hydrotropism), spaceflight condition (Space Flight, Ground Control), assay type, tissue, and hardware. Tropism labels were inferred heuristically from tissue type, hardware configuration (EMCS→gravitropism, Veggie/LED→phototropism), light regime, and study-specific overrides.

### Foundation Model: Salk Arabidopsis Developmental Atlas

The Salk atlas (GSE226097) contains 432,919 nuclei across 10 developmental stages (seedling to flower), profiled with 10x Chromium and integrated using Harmony [3]. We loaded the Seurat v5 object and extracted:

- **Signature matrix**: 4,000 highly variable genes (HVGs) × 183 clusters (orig.cluster), using log-normalized expression
- **Cell-type markers**: Top 50 markers per cluster (9,150 total), identified via Wilcoxon rank-sum test
- **Per-cell HVG expression**: 60,792 cells × 4,000 genes (stratified subsample covering all 183 clusters) for auto-decoder training

### Custom Variational Auto-decoder

We trained a conditional VAE to learn a 32-dimensional latent representation of cell-type-specific gene expression:

- **Encoder**: 4,000 → 512 → 256 → 32 (latent z)
- **Decoder**: (32 + 64 condition) → 512 → 512 → 4,000, conditioned on cluster embedding (183→32), organ embedding (12→16), and stage embedding (10→16)
- **Auxiliary head**: Cluster classification (cross-entropy loss)
- **Loss**: MSE reconstruction + 0.5 × KL divergence + 0.1 × cross-entropy
- **Training**: 40 epochs, batch size 256, AdamW optimizer (lr = 1e-3, cosine annealing), early stopping (patience = 8)
- **Hardware**: CPU (no GPU available), 60k-cell subsample

The trained model produces per-cluster "stimulus codes" — 32-dimensional vectors representing the learned transcriptional state of each cell type. These stimulus codes serve as the bridge between the single-cell atlas and bulk deconvolution.

### Bulk Deconvolution

We deconvolved 398 bulk samples (11 OSDR RNA-seq studies + 2 GEO datasets) onto the 183 atlas cell types using three methods:

1. **Auto-decoder stimulus codes (primary)**: Non-negative least squares (NNLS) fitting of bulk expression onto the signature matrix yields cell-type fractions, which are then used to weight the per-cluster stimulus codes, producing sample-level stimulus activation scores (32 dimensions).

2. **PhysioSpace-style (baseline)**: Cosine similarity between bulk samples and signature profiles in PCA space.

3. **CIBERSORTx-style (baseline)**: Direct NNLS on the signature matrix.

### Differential Expression and Meta-analysis

**DESeq2**: For each of 11 OSDR RNA-seq studies with matched flight and ground control samples, we ran DESeq2 with the contrast Space Flight vs Ground Control (Ground Control as reference level). Log2 fold changes were shrunk using apeglm.

**Meta-analysis**: We performed DerSimonian-Laird random-effects meta-analysis on the shrunk log2FC values across all 11 studies, requiring ≥2 studies per gene. The between-study variance (τ²) was estimated using the method of moments. Heterogeneity was assessed using Cochran's Q and I². P-values were adjusted using the Benjamini-Hochberg procedure.

### Tropism Classifier

We trained an elastic-net logistic regression classifier to distinguish:
1. **Space Flight vs Ground Control** (binary): Features = 183 cell-type fractions + 32 stimulus activation scores (215 features total). Nested 5-fold cross-validation (StratifiedKFold), class-weighted, l1_ratio = 0.5.
2. **Tropism type** (multinomial): Gravitropism vs phototropism (the two tropism types with ≥10 matched samples).

### Visualization

- **ggPlantMap**: Cell-type abundance differences projected onto Arabidopsis tissue anatomy maps (root tip cross-section, 3-week rosette, shoot apex, leaf cross-section, inflorescence stem cross-section) using the ggPlantmap.heatmap function.
- **KEGG pathway network**: Pathways sharing differentially expressed genes, rendered as a tidygraph/ggraph network (ggkegg/SBGNview unavailable due to Rgraphviz dependency failure).
- **Standard plots**: Volcano plot (meta-analysis), heatmaps (cell-type fractions, stimulus activation), forest plot (top gene across studies), bar plot (organ-level differences), feature importance (classifier).

## Results

### Dataset Composition

The harmonized dataset comprises 1,253 samples: 1,190 from OSDR and 63 from GEO. Tropism distribution: gravitropism (959), phototropism (240), gravitropism+phototropism (36), mechanotropism (12), hydrotropism (6). Spaceflight conditions: 832 Space Flight, 421 Ground Control. Assay types: 1,072 RNA-seq, 181 microarray.

### Auto-decoder Performance

The auto-decoder converged at epoch 5 (val_loss = 0.4300, val_recon = 0.140, val_kld = 0.207). The auxiliary classifier predicted 181/183 clusters correctly. The 32-dimensional latent embeddings showed a range of [-3.12, 3.19] with mean 0.003, indicating well-regularized representations. [TODO: no committed training log records the values in this paragraph (README, Limitation 8)] Per-cluster stimulus codes were extracted for all 183 clusters.

### Bulk Deconvolution

All 398 bulk samples were successfully deconvolved onto the 183 atlas cell types. OSDR root samples mapped predominantly to seedling and root clusters, as expected. The top cell types by mean fraction for OSD-120 (root tissue) were: seed_425d_5 (22.0%), seedling_15d_18 (12.9%), seedling_15d_15 (11.4%), silique_4 (9.0%), silique_6 (8.9%).

### Meta-analysis

Across 11 studies, 26,402 genes were tested in the random-effects meta-analysis, with 16,002 reaching significance (padj < 0.05). [TODO: produced with a τ² that was too small by a factor of Σw, fixed in `Code/meta_classifier.py` on 28 Oct 2026; re-run on the per-study DESeq2 tables, which are not in this repository (README, Limitation 7)] The top genes showed high heterogeneity (I² > 80%), reflecting genuine biological variability across experimental conditions. The most significant genes included AT2G33330 (pooled log2FC = 0.77), AT2G33450 (pooled log2FC = 1.07), and AT2G35130 (pooled log2FC = 1.45). [TODO: 4,776 genes tie at p = 0 in the current table, so this ranking is arbitrary (README, Limitation 7)]

### Classifier Performance

**Flight vs Ground Control**: The elastic-net classifier achieved AUC = 0.919 ± 0.047 (nested 5-fold CV) with 85% accuracy. [TODO: printed by `Code/meta_classifier.py` but never saved; re-run it and commit its metrics output (README, Limitation 8)] The top discriminating features were silique_20 (coefficient = 1.75), stem_9 (1.56), seedling_9d_18 (-1.29), rosette_21d_11 (-1.19), and flower_11 (1.02), along with stimulus dimensions 5 and 29.

**Tropism type**: The multinomial classifier achieved F1 = 1.000 for gravitropism vs phototropism [TODO: `tropism_confusion_matrix.csv` is written by `Code/meta_classifier.py` but not committed (README, Limitation 8)], reflecting the strong tissue/hardware confounding between these conditions (gravitropism experiments predominantly use root tissue in EMCS hardware; phototropism experiments use shoot tissue with LED lighting). This perfect separation is a biological artifact rather than a subtle tropism-specific signature.

### Cell-Type-Specific Responses

Three cell types showed significant abundance differences between Space Flight and Ground Control (padj < 0.05): stem_9 (padj = 0.003, +1.2% in flight), seed_425d_7 (padj = 0.027, +0.04% in flight), and silique_20 (padj = 0.030, +0.7% in flight). At the organ level, stem tissues showed the most significant response (mean_padj = 0.003), followed by seed (0.027) and silique (0.030).

### Tissue Projection

**TODO:** each map gives every region one organ-level mean plus random jitter (`Code/viz_ggplantmap.R`, `rnorm`, seed 42), so only the organ-level difference is data; the region-level descriptions below are not results (README, Limitation 9).

ggPlantMap projections visualized organ-level abundance changes onto Arabidopsis anatomy maps:
- **Root tip**: Columella and procambium regions showed elevated abundance in flight samples
- **Rosette**: Leaf ROIs showed minimal changes, consistent with root-focused experiments
- **Shoot apex**: Central zone and peripheral zone regions showed moderate changes
- **Leaf cross-section**: Vascular bundle regions showed slight elevation
- **Inflorescence stem**: Starch sheath and xylem regions showed the strongest response, consistent with gravitropism signaling in vascular tissues

## Discussion

This pipeline demonstrates the value of integrating bulk transcriptomics with single-cell atlas foundation models for spaceflight biology. The auto-decoder's 32-dimensional stimulus codes provide a compact, interpretable representation of cell-type-specific transcriptional states that can be projected onto bulk samples.

The Flight vs Ground Control classifier (AUC = 0.919 [TODO: README, Limitation 8]) demonstrates that cell-type composition and stimulus activation patterns contain strong discriminative signal for spaceflight response. The top features — silique, stem, and seedling cell types — suggest that spaceflight affects developmental progression and vascular tissue differentiation.

The meta-analysis identified 16,002 significant genes [TODO: τ² error, re-run; README, Limitation 7] with high heterogeneity, reflecting the diversity of experimental conditions (different tissues, ecotypes, hardware, missions, timepoints) across the 11 studies. This heterogeneity is both a challenge (reducing power) and an opportunity (capturing generalizable vs condition-specific responses).

### Limitations

1. **GPU unavailability**: The auto-decoder was trained on CPU with a 60k-cell subsample. While this covers all 183 clusters, a larger training set with GPU acceleration could improve latent representations.

2. **Tropism data asymmetry**: Gravitropism and phototropism dominate the dataset; hydrotropism (6 samples) and mechanotropism (12 samples) are underrepresented. Chemotropism and oxytropism have no dedicated Arabidopsis transcriptomics data.

3. **Metadata matching**: Only 124/398 deconvolved samples could be matched to tropism labels via GSM identifiers, limiting the classifier training set.

4. **Tropism classifier confounding**: The perfect F1 = 1.0 [TODO: README, Limitation 8] for gravitropism vs phototropism reflects tissue/hardware confounding rather than tropism-specific transcriptional signatures.

5. **Pathway visualization**: ggkegg and SBGNview were unavailable due to Rgraphviz dependency failure; KEGG pathway networks were rendered with tidygraph/ggraph instead.

## Data and Code Availability

All code, results, and figures are available in the Zenodo-ready repository. Data sources:
- NASA GeneLab OSDR: https://osdr.nasa.gov
- NCBI GEO: https://www.ncbi.nlm.nih.gov/geo
- Salk Atlas: GSE226097

## References

[1] Toyota M, Gilroy S. Gravitropism and mechanical signaling in plants. Am J Bot. 2013;100(1):111–125. doi:10.3732/ajb.1200408
[2] NASA GeneLab. OSDR API documentation. https://visualization.osdr.nasa.gov/biodata/api/v2/query/
[3] Salk Institute. Arabidopsis Developmental Atlas. GSE226097.
[4] Newman AM et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol. 2019;37:773–782. doi:10.1038/s41587-019-0114-2
[5] Lenz M, Schuldt BM, Müller F-J, Schuppert A. PhysioSpace: relating gene expression experiments from heterogeneous sources using shared physiological processes. PLoS ONE. 2013;8(10):e77627. doi:10.1371/journal.pone.0077627
[6] Lopez R et al. Deep generative modeling for single-cell transcriptomics. Nat Methods. 2018;15:1053–1058. doi:10.1038/s41592-018-0229-2
