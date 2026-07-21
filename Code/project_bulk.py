#!/usr/bin/env python3
"""
Subsystem 3 (continued): Bulk projection / cell-type deconvolution

Projects OSDR + GEO bulk RNA-seq samples onto the Salk atlas cell types using:
1. Auto-decoder stimulus codes (primary method) - non-negative least squares (NNLS)
   fitting of bulk expression onto the learned per-cluster stimulus code matrix
2. CIBERSORTx-style NNLS on signature matrix (baseline)
3. PhysioSpace-style projection (baseline) - cosine similarity in PCA space

Outputs per-sample cell-type fractions and stimulus activation scores.
"""
import os, json, gzip, glob
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

# Paths
ATLAS_DIR = '/mnt/shared-workspace/autodecoder'
PROC_DIR = '/mnt/shared-workspace/processed'
DE_DIR = f'{PROC_DIR}/de_results'
OUT_DIR = f'{PROC_DIR}/projection'
os.makedirs(OUT_DIR, exist_ok=True)

# ---- Load auto-decoder artifacts ----
print("=== Loading auto-decoder artifacts ===")
with open(f'{ATLAS_DIR}/model_encodings.json') as f:
    enc = json.load(f)
with open(f'{ATLAS_DIR}/cluster_stimulus_codes.json') as f:
    stim_codes = json.load(f)

# Build stimulus code matrix (n_clusters x latent_dim)
clusters = sorted(stim_codes.keys(), key=lambda x: enc['cluster_to_idx'].get(x, 999))
stim_matrix = np.array([stim_codes[c] for c in clusters])
print(f"Stimulus matrix: {stim_matrix.shape} (clusters x latent_dim)")

# Load HVGs
with open(f'{ATLAS_DIR}/hvgs.txt') as f:
    hvgs = [l.strip() for l in f if l.strip()]
print(f"HVGs: {len(hvgs)}")

# Load signature matrix (baseline method)
print("\n=== Loading signature matrix ===")
sig = pd.read_csv(f'{PROC_DIR}/cell_type_signatures.csv', index_col=0)
print(f"Signature matrix: {sig.shape} (genes x clusters)")
sig_clusters = list(sig.columns)
sig_genes = list(sig.index)

# Load cluster summary for organ/stage annotations
cluster_summary = pd.read_csv(f'{PROC_DIR}/cluster_summary.csv')
print(f"Cluster summary: {cluster_summary.shape}")
print(cluster_summary.columns.tolist())

# ---- Load bulk expression data ----
print("\n=== Loading bulk expression data ===")

def load_normalized_counts(study_id):
    """Load DESeq2-normalized counts for a study."""
    fpath = f'{DE_DIR}/{study_id}_normalized_counts.csv'
    if not os.path.exists(fpath):
        return None
    df = pd.read_csv(fpath, index_col=0)
    return df

def load_geo_counts():
    """Load GEO count matrices."""
    geo_data = {}
    # GSE143760 (phototropism) - featureCounts
    fpath = '/mnt/shared-workspace/raw_geo/GSE143760/GSE143760_counts.csv'
    if os.path.exists(fpath):
        df = pd.read_csv(fpath, index_col=0)
        df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
        lib_sizes = df.sum(axis=0).replace(0, 1)
        cpm = df.div(lib_sizes, axis=1) * 1e6
        geo_data['GSE143760'] = np.log1p(cpm).replace([np.inf, -np.inf], 0).fillna(0)
        print(f"  GSE143760: {df.shape}")
    # GSE225299 (mechanotropism) - pseudo-bulk from scRNA-seq
    fpath = '/mnt/shared-workspace/raw_geo/GSE225299/GSE225299_pseudobulk_counts.csv'
    if os.path.exists(fpath):
        df = pd.read_csv(fpath)
        if 'gene' in df.columns:
            df = df.set_index('gene')
        df = df.apply(pd.to_numeric, errors='coerce').fillna(0)
        lib_sizes = df.sum(axis=0).replace(0, 1)
        cpm = df.div(lib_sizes, axis=1) * 1e6
        geo_data['GSE225299'] = np.log1p(cpm).replace([np.inf, -np.inf], 0).fillna(0)
        print(f"  GSE225299: {df.shape}")
    return geo_data

# Collect all bulk data
bulk_datasets = {}
# OSDR studies with normalized counts
for f in glob.glob(f'{DE_DIR}/*_normalized_counts.csv'):
    sid = os.path.basename(f).replace('_normalized_counts.csv', '')
    df = pd.read_csv(f)
    # The gene names weren't written as index by R fwrite - get them from DE results
    de_file = f'{DE_DIR}/{sid}_de_results.tsv'
    if os.path.exists(de_file):
        de = pd.read_csv(de_file, sep='\t')
        if len(de) == len(df):
            df.index = de['gene'].values
        else:
            # Fallback: use AT gene pattern from first column if present
            print(f"  WARNING: {sid} row count mismatch ({len(de)} DE vs {len(df)} NC), using generic index")
            df.index = [f'gene_{i}' for i in range(len(df))]
    else:
        df.index = [f'gene_{i}' for i in range(len(df))]
    df = df.apply(pd.to_numeric, errors='coerce').dropna(how='all')
    bulk_datasets[sid] = np.log1p(df)
    print(f"  {sid}: {df.shape}")

# GEO data
geo_data = load_geo_counts()
bulk_datasets.update(geo_data)

print(f"\nTotal bulk datasets: {len(bulk_datasets)}")
total_samples = sum(df.shape[1] for df in bulk_datasets.values())
print(f"Total bulk samples: {total_samples}")

# ---- Deconvolution methods ----

def deconvolve_nnls_signature(bulk_df, sig_df, min_overlap=1000):
    """CIBERSORTx-style: NNLS fitting of bulk onto signature matrix.
    Returns per-sample cell-type fractions."""
    common_genes = list(set(bulk_df.index) & set(sig_df.index))
    if len(common_genes) < min_overlap:
        print(f"  WARNING: only {len(common_genes)} overlapping genes")
        if len(common_genes) < 100:
            return None

    B = bulk_df.loc[common_genes].values  # genes x samples
    S = sig_df.loc[common_genes].values   # genes x clusters

    # Replace any inf/nan with 0
    B = np.nan_to_num(B, nan=0.0, posinf=0.0, neginf=0.0)
    S = np.nan_to_num(S, nan=0.0, posinf=0.0, neginf=0.0)

    fractions = np.zeros((bulk_df.shape[1], S.shape[1]))
    for j in range(bulk_df.shape[1]):
        # NNLS: solve B[:,j] ≈ S @ w, w >= 0
        w, _ = nnls(S, B[:, j])
        # Normalize to fractions
        if w.sum() > 0:
            fractions[j, :] = w / w.sum()

    return pd.DataFrame(fractions, index=bulk_df.columns, columns=sig_df.columns)

def deconvolve_stimulus_codes(bulk_df, sig_df, stim_matrix, clusters, min_overlap=1000):
    """Auto-decoder method: Project bulk onto signature space, then map to stimulus codes.
    Step 1: NNLS to get cell-type proportions from signature matrix
    Step 2: Weight stimulus codes by proportions to get sample-level stimulus activation
    """
    # First get proportions via NNLS on signature
    frac_df = deconvolve_nnls_signature(bulk_df, sig_df, min_overlap)
    if frac_df is None:
        return None, None

    # Map signature clusters to stimulus code clusters
    # stim_matrix is ordered by `clusters` list
    cluster_to_stim_idx = {c: i for i, c in enumerate(clusters)}

    # Compute stimulus activation: weighted average of stimulus codes
    # frac_df columns are signature clusters; need to align
    common_clusters = [c for c in frac_df.columns if c in cluster_to_stim_idx]
    if len(common_clusters) < len(frac_df.columns) * 0.5:
        print(f"  WARNING: only {len(common_clusters)}/{len(frac_df.columns)} clusters matched to stimulus codes")
        # Try fuzzy matching
        common_clusters = []
        for sc in frac_df.columns:
            # Try exact, then prefix match
            if sc in cluster_to_stim_idx:
                common_clusters.append(sc)
            else:
                # Try matching by stripping suffixes
                for tc in cluster_to_stim_idx:
                    if tc.startswith(sc) or sc.startswith(tc):
                        common_clusters.append(tc)
                        break

    if len(common_clusters) == 0:
        print("  ERROR: no cluster overlap between signature and stimulus codes")
        return frac_df, None

    # Stimulus activation matrix: samples x latent_dim
    stim_indices = [cluster_to_stim_idx[c] for c in common_clusters]
    frac_subset = frac_df[common_clusters].values
    stim_activation = frac_subset @ stim_matrix[stim_indices, :]

    stim_df = pd.DataFrame(stim_activation,
                           index=frac_df.index,
                           columns=[f'stim_dim_{i}' for i in range(stim_matrix.shape[1])])
    return frac_df, stim_df

def physiospace_projection(bulk_df, sig_df, min_overlap=1000):
    """PhysioSpace-style: cosine similarity in PCA space.
    Projects bulk samples onto reference (signature) and computes similarity scores."""
    common_genes = list(set(bulk_df.index) & set(sig_df.index))
    if len(common_genes) < min_overlap:
        return None

    B = bulk_df.loc[common_genes].values  # genes x samples
    S = sig_df.loc[common_genes].values   # genes x clusters

    # Replace any inf/nan with 0
    B = np.nan_to_num(B, nan=0.0, posinf=0.0, neginf=0.0)
    S = np.nan_to_num(S, nan=0.0, posinf=0.0, neginf=0.0)

    # Combine and PCA
    combined = np.hstack([S, B])
    pca = PCA(n_components=min(50, min(combined.shape)))
    combined_pca = pca.fit_transform(combined.T)  # (n_clusters + n_samples) x PCs

    S_pca = combined_pca[:S.shape[1], :]  # clusters x PCs
    B_pca = combined_pca[S.shape[1]:, :]  # samples x PCs

    # Cosine similarity
    from sklearn.metrics.pairwise import cosine_similarity
    sim = cosine_similarity(B_pca, S_pca)  # samples x clusters

    # Normalize to pseudo-fractions (softmax-like)
    sim_norm = (sim - sim.min(axis=1, keepdims=True)) / (sim.max(axis=1, keepdims=True) - sim.min(axis=1, keepdims=True) + 1e-8)
    return pd.DataFrame(sim_norm, index=bulk_df.columns, columns=sig_df.columns)

# ---- Run deconvolution on all datasets ----
print("\n=== Running deconvolution ===")
all_fractions = {}
all_stimulus = {}
all_physiospace = {}

for ds_name, bulk_df in bulk_datasets.items():
    print(f"\n--- {ds_name} ({bulk_df.shape[1]} samples) ---")

    # Ensure gene IDs are strings
    bulk_df.index = bulk_df.index.astype(str)
    sig.index = sig.index.astype(str)

    # Method 1: Auto-decoder stimulus codes
    frac_df, stim_df = deconvolve_stimulus_codes(bulk_df, sig, stim_matrix, clusters)
    if frac_df is not None:
        all_fractions[ds_name] = frac_df
        print(f"  NNLS fractions: {frac_df.shape}")
        if stim_df is not None:
            all_stimulus[ds_name] = stim_df
            print(f"  Stimulus activation: {stim_df.shape}")

    # Method 2: PhysioSpace baseline
    phys_df = physiospace_projection(bulk_df, sig)
    if phys_df is not None:
        all_physiospace[ds_name] = phys_df
        print(f"  PhysioSpace: {phys_df.shape}")

# ---- Save results ----
print("\n=== Saving results ===")

# Combine all fractions
frac_combined = []
for ds_name, frac_df in all_fractions.items():
    frac_df_copy = frac_df.copy()
    frac_df_copy['dataset'] = ds_name
    frac_df_copy['sample'] = frac_df.index
    frac_combined.append(frac_df_copy)
if frac_combined:
    frac_all = pd.concat(frac_combined, ignore_index=True)
    frac_all.to_csv(f'{OUT_DIR}/cell_type_fractions_all.csv', index=False)
    print(f"Saved fractions: {frac_all.shape}")

# Combine stimulus activations
stim_combined = []
for ds_name, stim_df in all_stimulus.items():
    stim_df_copy = stim_df.copy()
    stim_df_copy['dataset'] = ds_name
    stim_df_copy['sample'] = stim_df.index
    stim_combined.append(stim_df_copy)
if stim_combined:
    stim_all = pd.concat(stim_combined, ignore_index=True)
    stim_all.to_csv(f'{OUT_DIR}/stimulus_activation_all.csv', index=False)
    print(f"Saved stimulus activations: {stim_all.shape}")

# Combine PhysioSpace
phys_combined = []
for ds_name, phys_df in all_physiospace.items():
    phys_df_copy = phys_df.copy()
    phys_df_copy['dataset'] = ds_name
    phys_df_copy['sample'] = phys_df.index
    phys_combined.append(phys_df_copy)
if phys_combined:
    phys_all = pd.concat(phys_combined, ignore_index=True)
    phys_all.to_csv(f'{OUT_DIR}/physiospace_scores_all.csv', index=False)
    print(f"Saved PhysioSpace scores: {phys_all.shape}")

# ---- Summary statistics ----
print("\n=== Summary ===")
print(f"Datasets deconvolved: {len(all_fractions)}")
print(f"Total samples with fractions: {sum(df.shape[0] for df in all_fractions.values())}")
print(f"Total samples with stimulus codes: {sum(df.shape[0] for df in all_stimulus.values())}")

# Show example
if all_fractions:
    example_ds = list(all_fractions.keys())[0]
    example_frac = all_fractions[example_ds]
    print(f"\nExample ({example_ds}) - top 5 cell types by mean fraction:")
    mean_frac = example_frac.mean(axis=0).sort_values(ascending=False)
    print(mean_frac.head(5))

print("\nDONE_PROJECTION")
