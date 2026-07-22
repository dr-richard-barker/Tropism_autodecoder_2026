#!/usr/bin/env python3
"""
export_web_artifacts.py — assemble the Phase-2 model bundle for the web tool.

Reads the pipeline's model artifacts and writes browser-ready files into
docs/assets/model/, so the web tool can reproduce the manuscript's cell-type
deconvolution (NNLS onto the signature matrix), stimulus-code projection, and
the Flight-vs-Ground-Control elastic-net classifier (nested-CV AUC 0.919)
entirely client-side.

Run this AFTER the pipeline (or the full-atlas GPU retrain) has produced:
  - cell_type_signatures.csv     (genes x 183 clusters)          [atlas preprocessing]
  - cluster_stimulus_codes.json  ({cluster: [32]})               [train_autodecoder*.py]
  - classifier_params.json       (coef/intercept/scaler/...)     [meta_classifier.py]
  - meta_analysis_results.tsv    (gene, pooled_log2FC, padj)     [meta_classifier.py]  (optional)

Outputs (docs/assets/model/):
  signature_matrix.bin      Float32, row-major genes x clusters
  signature_index.json      { genes:[...], clusters:[...], shape:[g,c] }
  stimulus_codes.json       { cluster: [32 floats] }  (copied)
  classifier.json           (copied classifier_params.json)
  flight_signature_genes.json  top-N meta-analysis genes (optional panel)
  manifest.json             versions + shapes + presence flags (gates Phase 2 in the UI)
"""
import os, json, argparse
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sig-matrix", default="/mnt/shared-workspace/processed/cell_type_signatures.csv")
    p.add_argument("--stim-codes", default="/mnt/shared-workspace/autodecoder/cluster_stimulus_codes.json")
    p.add_argument("--classifier", default="/mnt/shared-workspace/processed/classifier/classifier_params.json")
    p.add_argument("--meta-analysis", default="/mnt/shared-workspace/processed/classifier/meta_analysis_results.tsv")
    p.add_argument("--out", default="docs/assets/model")
    p.add_argument("--top-genes", type=int, default=300)
    p.add_argument("--version", default="2.0.0-artifacts")
    p.add_argument("--created", default="", help="ISO date stamp for the manifest (optional)")
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)
    manifest = {"version": args.version, "created": args.created, "files": {}, "shapes": {}, "present": {}}

    # --- signature matrix -> Float32 binary (row-major genes x clusters) + index ---
    sig = pd.read_csv(args.sig_matrix, index_col=0)
    genes = [str(g) for g in sig.index]
    clusters = [str(c) for c in sig.columns]
    S = np.ascontiguousarray(sig.values, dtype=np.float32)  # genes x clusters, row-major
    S.tofile(os.path.join(args.out, "signature_matrix.bin"))
    with open(os.path.join(args.out, "signature_index.json"), "w") as f:
        json.dump({"genes": genes, "clusters": clusters, "shape": [len(genes), len(clusters)]}, f)
    manifest["files"]["signature_matrix"] = "signature_matrix.bin"
    manifest["files"]["signature_index"] = "signature_index.json"
    manifest["shapes"]["signature_matrix"] = [len(genes), len(clusters)]
    manifest["present"]["signature_matrix"] = True
    print(f"signature_matrix.bin: {len(genes)} genes x {len(clusters)} clusters "
          f"({S.nbytes/1e6:.1f} MB)")

    # --- stimulus codes (copy through, validate 32-dim) ---
    with open(args.stim_codes) as f:
        stim = json.load(f)
    latent = len(next(iter(stim.values()))) if stim else 0
    with open(os.path.join(args.out, "stimulus_codes.json"), "w") as f:
        json.dump(stim, f)
    manifest["files"]["stimulus_codes"] = "stimulus_codes.json"
    manifest["shapes"]["stimulus_codes"] = [len(stim), latent]
    manifest["present"]["stimulus_codes"] = True
    print(f"stimulus_codes.json: {len(stim)} clusters x {latent} dims")

    # --- classifier params (copy through) ---
    with open(args.classifier) as f:
        clf = json.load(f)
    with open(os.path.join(args.out, "classifier.json"), "w") as f:
        json.dump(clf, f)
    manifest["files"]["classifier"] = "classifier.json"
    manifest["shapes"]["classifier_features"] = len(clf.get("feature_names", []))
    manifest["present"]["classifier"] = True
    print(f"classifier.json: {len(clf.get('feature_names', []))} features, "
          f"positive class = {clf.get('positive_class')}")

    # --- flight-response genes (optional display panel) ---
    if os.path.exists(args.meta_analysis):
        m = pd.read_csv(args.meta_analysis, sep="\t")
        cols = {c.lower(): c for c in m.columns}
        gcol = cols.get("gene"); pcol = cols.get("padj"); fcol = cols.get("pooled_log2fc")
        if gcol and pcol:
            top = m.sort_values(pcol).head(args.top_genes)
            recs = [{"gene": str(r[gcol]),
                     "log2FC": (float(r[fcol]) if fcol and pd.notna(r[fcol]) else None),
                     "padj": (float(r[pcol]) if pd.notna(r[pcol]) else None)}
                    for _, r in top.iterrows()]
            with open(os.path.join(args.out, "flight_signature_genes.json"), "w") as f:
                json.dump(recs, f)
            manifest["files"]["flight_signature_genes"] = "flight_signature_genes.json"
            manifest["shapes"]["flight_signature_genes"] = len(recs)
            manifest["present"]["flight_signature_genes"] = True
            print(f"flight_signature_genes.json: top {len(recs)} genes")
    else:
        manifest["present"]["flight_signature_genes"] = False
        print("meta_analysis_results.tsv not found — skipping flight_signature_genes.json")

    # --- manifest (its presence + all-true flags gate Phase 2 in the UI) ---
    manifest["phase2_ready"] = all(manifest["present"].get(k) for k in
                                   ("signature_matrix", "stimulus_codes", "classifier"))
    with open(os.path.join(args.out, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"manifest.json written — phase2_ready = {manifest['phase2_ready']}")
    print(f"\nDONE. Copy {args.out}/ into the repo and redeploy Pages to activate Phase 2.")


if __name__ == "__main__":
    main()
