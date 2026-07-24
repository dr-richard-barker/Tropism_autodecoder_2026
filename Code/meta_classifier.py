#!/usr/bin/env python3
"""
Subsystem 4 (continued): Meta-analysis + Tropism classifier

1. Meta-analysis: Random-effects meta-analysis (DerSimonian-Laird) on log2FC
   across all OSDR studies, stratified by tropism type.
2. Tropism classifier: Elastic-net logistic regression with nested 5-fold CV
   to classify tropism type from cell-type fractions + stimulus activation.
3. Cell-type-specific tropism scoring: Identify cell types whose abundance
   or stimulus activation differs between Flight and Ground Control.
"""
import os, json, glob
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression, ElasticNet
from sklearn.model_selection import StratifiedKFold, cross_val_score, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.utils.class_weight import compute_class_weight
import warnings
warnings.filterwarnings('ignore')

import argparse
_ap = argparse.ArgumentParser(description="Meta-analysis + Flight-vs-GC / tropism classifier.")
_ap.add_argument('--meta', default=os.environ.get('META', 'Data/harmonized_metadata.tsv'))
_ap.add_argument('--de-dir', default=os.environ.get('DE_DIR', 'bulk/de_results'))
_ap.add_argument('--proj-dir', default=os.environ.get('PROJ_DIR', 'bulk/projection'))
_ap.add_argument('--out-dir', default=os.environ.get('OUT_DIR', 'bulk/classifier'))
_args, _ = _ap.parse_known_args()
META = _args.meta
DE_DIR = _args.de_dir
PROJ_DIR = _args.proj_dir
OUT_DIR = _args.out_dir
os.makedirs(OUT_DIR, exist_ok=True)

# ---- Load harmonized metadata ----
print("=== Loading metadata ===")
meta = pd.read_csv(META, sep='\t')
print(f"Metadata: {meta.shape}")
print(f"Tropism distribution:\n{meta['tropism_type'].value_counts()}")
print(f"\nCondition distribution:\n{meta['spaceflight_condition'].value_counts()}")

# ---- Load DE results ----
print("\n=== Loading DE results ===")
de_files = glob.glob(f'{DE_DIR}/*_de_results.tsv')
de_results = {}
for f in de_files:
    sid = os.path.basename(f).replace('_de_results.tsv', '')
    de = pd.read_csv(f, sep='\t')
    de_results[sid] = de
    print(f"  {sid}: {len(de)} genes, {de['padj'].notna().sum()} tested")

# ---- Meta-analysis (DerSimonian-Laird random effects) ----
print("\n=== Meta-analysis ===")

def random_effects_meta(effect_sizes, variances):
    """DerSimonian-Laird random-effects meta-analysis.
    Returns pooled effect, SE, z, p, I2, Q."""
    effect_sizes = np.array(effect_sizes, dtype=float)
    variances = np.array(variances, dtype=float)

    # Remove NaN/inf
    mask = np.isfinite(effect_sizes) & np.isfinite(variances) & (variances > 0)
    effect_sizes = effect_sizes[mask]
    variances = variances[mask]

    if len(effect_sizes) < 2:
        return None

    # Fixed-effect weights
    w_fe = 1.0 / variances
    mu_fe = np.sum(w_fe * effect_sizes) / np.sum(w_fe)

    # Cochran's Q
    Q = np.sum(w_fe * (effect_sizes - mu_fe) ** 2)
    df = len(effect_sizes) - 1

    # Between-study variance (tau^2)
    if Q > df:
        tau2 = (Q - df) / np.sum(w_fe) * (np.sum(w_fe) / (np.sum(w_fe) ** 2 - np.sum(w_fe ** 2)))
    else:
        tau2 = 0.0

    # Random-effects weights
    w_re = 1.0 / (variances + tau2)
    mu_re = np.sum(w_re * effect_sizes) / np.sum(w_re)
    se_re = np.sqrt(1.0 / np.sum(w_re))
    z = mu_re / se_re
    p = 2 * (1 - stats.norm.cdf(abs(z)))

    # I^2
    I2 = max(0, (Q - df) / Q * 100) if Q > 0 else 0.0

    return {
        'pooled_log2FC': mu_re, 'se': se_re, 'z': z, 'p': p,
        'tau2': tau2, 'Q': Q, 'I2': I2, 'n_studies': len(effect_sizes)
    }

# Get all genes across studies
all_genes = set()
for de in de_results.values():
    all_genes.update(de['gene'].tolist())
all_genes = sorted(all_genes)
print(f"Total unique genes across studies: {len(all_genes)}")

# Run meta-analysis per gene
meta_results = []
for gene in all_genes:
    l2fcs = []
    vars_ = []
    studies = []
    for sid, de in de_results.items():
        row = de[de['gene'] == gene]
        if len(row) > 0 and pd.notna(row.iloc[0]['padj']):
            l2fc = row.iloc[0]['log2FC_shrunk'] if pd.notna(row.iloc[0].get('log2FC_shrunk')) else row.iloc[0]['log2FC']
            se = row.iloc[0]['lfcSE']
            if pd.notna(l2fc) and pd.notna(se) and se > 0:
                l2fcs.append(l2fc)
                vars_.append(se ** 2)
                studies.append(sid)

    if len(l2fcs) >= 2:
        result = random_effects_meta(l2fcs, vars_)
        if result:
            result['gene'] = gene
            result['studies'] = ';'.join(studies)
            meta_results.append(result)

meta_df = pd.DataFrame(meta_results)
# BH correction
from statsmodels.stats.multitest import multipletests
meta_df['padj'] = multipletests(meta_df['p'], method='fdr_bh')[1]
meta_df = meta_df.sort_values('p')
meta_df.to_csv(f'{OUT_DIR}/meta_analysis_results.tsv', sep='\t', index=False)
print(f"\nMeta-analysis: {len(meta_df)} genes tested across {len(de_results)} studies")
print(f"Significant (padj<0.05): {(meta_df['padj'] < 0.05).sum()}")
print(f"Top 10 genes:")
print(meta_df[['gene','pooled_log2FC','p','padj','I2','n_studies']].head(10).to_string())

# ---- Load projection results and build classifier features ----
print("\n=== Building classifier features ===")
fractions = pd.read_csv(f'{PROJ_DIR}/cell_type_fractions_all.csv')
stimulus = pd.read_csv(f'{PROJ_DIR}/stimulus_activation_all.csv')

print(f"Fractions: {fractions.shape}")
print(f"Stimulus: {stimulus.shape}")

# Merge with metadata to get tropism labels and conditions
# The sample column in fractions/stimulus needs to match metadata
print(f"\nFraction samples (first 5): {fractions['sample'].head().tolist()}")
print(f"Metadata columns: {list(meta.columns)}")

# Build sample key for matching
# OSDR samples have format like "Atha_Col-0_root_GC_Alight_Rep1_GSM2493759_Day13"
# Metadata has sample_id, GSM, OSD study etc.
# Let's try matching on GSM if present in sample name
import re
def extract_gsm(s):
    m = re.search(r'(GSM\d+)', str(s))
    return m.group(1) if m else None

fractions['GSM'] = fractions['sample'].apply(extract_gsm)
stimulus['GSM'] = stimulus['sample'].apply(extract_gsm)

# Check metadata for GSM column
gsm_cols = [c for c in meta.columns if 'gsm' in c.lower() or 'sample' in c.lower() or 'accession' in c.lower()]
print(f"Metadata GSM-like columns: {gsm_cols}")

# Try to match
if 'GSM' in meta.columns:
    meta_gsm = meta.set_index('GSM')
elif 'sample_id' in meta.columns:
    meta_gsm = meta.set_index('sample_id')
else:
    # Use the first column that has GSM-like values
    for c in meta.columns:
        if meta[c].astype(str).str.contains('GSM').any():
            meta_gsm = meta.set_index(c)
            print(f"Using {c} as sample key")
            break

# Match fractions to metadata
matched = 0
for idx, row in fractions.iterrows():
    gsm = row['GSM']
    if gsm and gsm in meta_gsm.index:
        matched += 1
print(f"\nMatched {matched}/{len(fractions)} samples to metadata via GSM")

# For unmatched, try matching on sample name directly
if matched < len(fractions) * 0.5:
    # Try matching on dataset + sample
    print("Trying alternative matching...")
    # Build a combined key
    fractions['match_key'] = fractions['dataset'] + '_' + fractions['sample']
    # Check if metadata has study + sample info
    print(f"Metadata sample columns: {gsm_cols}")
    # For OSDR, the sample names contain GSM IDs
    # Let's check what fraction of metadata has GSM
    if 'GSM' in meta.columns:
        print(f"Metadata GSM coverage: {meta['GSM'].notna().sum()}/{len(meta)}")

# Build feature matrix from samples that matched
# Use GSM matching
if 'GSM' in meta.columns:
    meta_indexed = meta.dropna(subset=['GSM']).set_index('GSM')
else:
    meta_indexed = meta_gsm

# Merge
frac_with_labels = fractions.dropna(subset=['GSM']).merge(
    meta_indexed[['tropism_type', 'spaceflight_condition', 'assay_type', 'tissue', 'genotype']],
    left_on='GSM', right_index=True, how='inner'
)
stim_with_labels = stimulus.dropna(subset=['GSM']).merge(
    meta_indexed[['tropism_type', 'spaceflight_condition', 'assay_type', 'tissue', 'genotype']],
    left_on='GSM', right_index=True, how='inner'
)

print(f"\nMatched fractions with labels: {len(frac_with_labels)}")
print(f"Tropism distribution in matched:\n{frac_with_labels['tropism_type'].value_counts()}")

# ---- Tropism classifier ----
print("\n=== Tropism classifier ===")

# Feature matrix: cell-type fractions + stimulus activation
# Merge fractions and stimulus on sample
features = frac_with_labels.merge(
    stim_with_labels[['sample'] + [c for c in stim_with_labels.columns if c.startswith('stim_dim_')]],
    on='sample', how='inner'
)
print(f"Combined features: {features.shape}")

# Get feature columns (exclude metadata)
feature_cols = [c for c in features.columns if c not in
                ['sample', 'dataset', 'GSM', 'tropism_type', 'spaceflight_condition', 'assay_type', 'tissue', 'genotype', 'match_key']]

X = features[feature_cols].values
y_tropism = features['tropism_type'].values
y_condition = features['spaceflight_condition'].values

print(f"X: {X.shape}, y_tropism: {len(y_tropism)}, y_condition: {len(y_condition)}")

# Classifier 1: Flight vs Ground Control (binary)
print("\n--- Flight vs Ground Control ---")
# Filter to only Flight/Ground Control samples
mask = pd.Series(y_condition).isin(['Space Flight', 'Ground Control']).values
X_cond = X[mask]
y_cond = y_condition[mask]
print(f"Samples: {len(y_cond)} (Flight={sum(y_cond=='Space Flight')}, GC={sum(y_cond=='Ground Control')})")

if len(y_cond) > 20 and len(np.unique(y_cond)) == 2:
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_cond)

    # Elastic-net logistic regression with nested CV
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    clf = LogisticRegression(
        penalty='elasticnet', solver='saga', max_iter=5000,
        l1_ratio=0.5, C=1.0, class_weight='balanced'
    )

    cv_scores = cross_val_score(clf, X_scaled, y_cond, cv=outer_cv, scoring='roc_auc')
    print(f"Nested CV AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

    y_pred = cross_val_predict(clf, X_scaled, y_cond, cv=outer_cv)
    print(f"\nClassification report (Flight vs GC):")
    print(classification_report(y_cond, y_pred))

    # Fit final model and get feature importance
    clf.fit(X_scaled, y_cond)
    importance = pd.DataFrame({
        'feature': feature_cols,
        'coefficient': clf.coef_[0]
    }).sort_values('coefficient', key=abs, ascending=False)
    importance.to_csv(f'{OUT_DIR}/flight_vs_gc_feature_importance.csv', index=False)
    print(f"\nTop 10 features for Flight vs GC:")
    print(importance.head(10).to_string())

    # Export the FULL fitted classifier for the web tool (Phase 2 artifact).
    # The feature-importance CSV above omits the intercept and the StandardScaler
    # parameters, so it alone cannot reproduce predictions — dump everything here.
    with open(f'{OUT_DIR}/classifier_params.json', 'w') as f:
        json.dump({
            'target': 'flight_vs_gc',
            'classes': clf.classes_.tolist(),
            'positive_class': str(clf.classes_[1]),   # coef_[0] is for classes_[1]
            'feature_names': feature_cols,             # 183 fractions (sig order) + stim_dim_0..31
            'coef': clf.coef_[0].tolist(),
            'intercept': float(clf.intercept_[0]),
            'scaler_mean': scaler.mean_.tolist(),
            'scaler_scale': scaler.scale_.tolist(),
            'l1_ratio': 0.5, 'C': 1.0,
            'n_train': int(len(y_cond)),
            'cv_auc_mean': float(cv_scores.mean()), 'cv_auc_std': float(cv_scores.std()),
        }, f, indent=2)
    print(f"Wrote classifier_params.json ({len(feature_cols)} features, "
          f"positive class = {clf.classes_[1]})")

# Classifier 2: Tropism type (multiclass)
print("\n--- Tropism type classifier ---")
# Filter to tropism types with enough samples
tropism_counts = pd.Series(y_tropism).value_counts()
valid_tropisms = tropism_counts[tropism_counts >= 10].index.tolist()
print(f"Tropism types with >=10 samples: {valid_tropisms}")

mask = pd.Series(y_tropism).isin(valid_tropisms).values
X_trop = X[mask]
y_trop = y_tropism[mask]
print(f"Samples: {len(y_trop)}")

if len(y_trop) > 30 and len(valid_tropisms) >= 2:
    scaler_t = StandardScaler()
    X_scaled_t = scaler_t.fit_transform(X_trop)

    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    clf_t = LogisticRegression(
        penalty='elasticnet', solver='saga', max_iter=5000,
        l1_ratio=0.5, C=1.0, class_weight='balanced',
        multi_class='multinomial'
    )

    cv_scores_t = cross_val_score(clf_t, X_scaled_t, y_trop, cv=outer_cv, scoring='f1_macro')
    print(f"Nested CV F1-macro: {cv_scores_t.mean():.3f} ± {cv_scores_t.std():.3f}")

    y_pred_t = cross_val_predict(clf_t, X_scaled_t, y_trop, cv=outer_cv)
    print(f"\nClassification report (tropism type):")
    print(classification_report(y_trop, y_pred_t))

    # Confusion matrix
    cm = confusion_matrix(y_trop, y_pred_t, labels=valid_tropisms)
    cm_df = pd.DataFrame(cm, index=valid_tropisms, columns=valid_tropisms)
    cm_df.to_csv(f'{OUT_DIR}/tropism_confusion_matrix.csv')

    # Fit final model
    clf_t.fit(X_scaled_t, y_trop)
    # Feature importance per tropism
    for i, trop in enumerate(clf_t.classes_):
        if i >= len(clf_t.coef_): break
        imp = pd.DataFrame({
            'feature': feature_cols,
            'coefficient': clf_t.coef_[i]
        }).sort_values('coefficient', key=abs, ascending=False)
        imp.to_csv(f'{OUT_DIR}/tropism_{trop}_feature_importance.csv', index=False)
        print(f"\nTop 5 features for {trop}:")
        print(imp.head(5).to_string())

# ---- Cell-type-specific tropism scoring ----
print("\n=== Cell-type-specific tropism scoring ===")

# Compare cell-type fractions between Flight and Ground Control
frac_labeled = frac_with_labels.copy()
frac_labeled = frac_labeled[frac_labeled['spaceflight_condition'].isin(['Space Flight', 'Ground Control'])]

cluster_cols = [c for c in frac_labeled.columns if c not in
                ['sample', 'dataset', 'GSM', 'tropism_type', 'spaceflight_condition', 'assay_type', 'tissue', 'genotype', 'match_key']]

flight = frac_labeled[frac_labeled['spaceflight_condition'] == 'Space Flight']
ground = frac_labeled[frac_labeled['spaceflight_condition'] == 'Ground Control']

print(f"Flight: {len(flight)}, Ground: {len(ground)}")

celltype_diff = []
for col in cluster_cols:
    f_vals = flight[col].dropna()
    g_vals = ground[col].dropna()
    if len(f_vals) > 5 and len(g_vals) > 5:
        stat, p = stats.mannwhitneyu(f_vals, g_vals, alternative='two-sided')
        mean_diff = f_vals.mean() - g_vals.mean()
        celltype_diff.append({
            'cell_type': col,
            'mean_flight': f_vals.mean(),
            'mean_ground': g_vals.mean(),
            'mean_diff': mean_diff,
            'pvalue': p,
            'n_flight': len(f_vals),
            'n_ground': len(g_vals)
        })

ct_diff_df = pd.DataFrame(celltype_diff)
if len(ct_diff_df) > 0:
    ct_diff_df['padj'] = multipletests(ct_diff_df['pvalue'], method='fdr_bh')[1]
    ct_diff_df = ct_diff_df.sort_values('padj')
    ct_diff_df.to_csv(f'{OUT_DIR}/celltype_flight_vs_ground.csv', index=False)
    print(f"\nCell types tested: {len(ct_diff_df)}")
    print(f"Significant (padj<0.05): {(ct_diff_df['padj'] < 0.05).sum()}")
    print(f"\nTop 10 differential cell types:")
    print(ct_diff_df[['cell_type','mean_flight','mean_ground','mean_diff','padj']].head(10).to_string())

print("\nDONE_CLASSIFIER")
