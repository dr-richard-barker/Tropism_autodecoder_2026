# Where each file in `Results/` comes from

Recorded 28 Oct 2026. "v1.0.0" means the July 2026 release; its intermediate workspace (per-study DESeq2 tables,
training log) was not retained.

| File | Produced by | Inputs | Notes |
| :- | :- | :- | :- |
| `cell_type_fractions_all.csv`, `stimulus_activation_all.csv`, `physiospace_scores_all.csv` | `Code/project_bulk.py` (v1.0.0) | v1.0.0 bulk counts and atlas signatures (not retained) | Unchanged since v1.0.0. |
| `celltype_flight_vs_ground.csv` | `Code/meta_classifier.py` | the two projection tables above, `Data/harmonized_metadata.tsv` | Re-run on 28 Oct 2026 reproduces the committed table to within 1e-16. |
| `classifier_metrics.json`, `classifier_params.json`, `tropism_confusion_matrix.csv`, `tropism_gravitropism_feature_importance.csv`, `flight_vs_gc_feature_importance.csv` | `Code/meta_classifier.py`, re-run 28 Oct 2026 | the two projection tables above, `Data/harmonized_metadata.tsv` | Reproduces the v1.0.0 numbers: AUC 0.919 ± 0.047, accuracy 0.85, tropism F1 1.000 (80 phototropism, 44 gravitropism). Coefficients agree with the v1.0.0 table to within 0.0008. Package versions are in `classifier_metrics.json`. |
| `meta_analysis_results.tsv` | `Code/meta_classifier.py` with the corrected DerSimonian–Laird τ² and p-value (28 Oct 2026) | DESeq2 results for the 11 studies in `de_analysis_summary.tsv`, from OSDR counts re-acquired in July 2026 | 26,527 genes, 248 with padj < 0.05. The v1.0.0 table (26,402 genes, 16,002 significant) came from the uncorrected code and the lost v1.0.0 DESeq2 tables. With the uncorrected code the re-acquired inputs give 15,926 significant genes. The re-acquired DESeq2 tables (20 MB) are not in git [TODO: deposit on Zenodo]. |
| `de_analysis_summary.tsv` | `Code/run_de.R` (the local copy differed only in a hard-coded path) | counts from `Code/fetch_osdr_counts.py` → `Code/build_de_inputs.py`, July 2026 | The 11 studies used by the meta-analysis. Gene and DE counts differ slightly from v1.0.0 because the counts were re-acquired. |
| `autodecoder_fullatlas_training_log.csv`, `autodecoder_fullatlas_hparams.yaml` | `Code/train_autodecoder_gpu.py` (Lightning `CSVLogger`), 23 Jul 2026 | full 432,919-nucleus atlas export (`Code/export_atlas.R`) | Full-atlas retrain: best val_loss 0.3739 (epoch 3), val_recon 0.139, val_kld 0.201. No log survives for the v1.0.0 draft model. |

To reproduce the classifier and meta-analysis files:

```bash
python Code/meta_classifier.py --meta Data/harmonized_metadata.tsv \
       --de-dir <dir with the 11 *_de_results.tsv> \
       --proj-dir <dir with cell_type_fractions_all.csv and stimulus_activation_all.csv> \
       --out-dir Results
```
