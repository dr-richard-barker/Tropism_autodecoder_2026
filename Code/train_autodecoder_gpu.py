#!/usr/bin/env python3
"""Subsystem 3 (GPU): Custom stimulus auto-decoder (conditional VAE) — full-atlas / GPU build.

This is a GPU-enabled, memory-scalable refactor of Code/train_autodecoder.py. The MODEL
(StimulusAutoDecoder) and LOSS are unchanged so results remain directly comparable to the
CPU/60k first draft; only the data path, device handling, and CLI change.

What changed vs the original
----------------------------
1. Device: auto-selects GPU (accelerator="gpu", configurable) with mixed precision
   (bf16-mixed on Ampere+; falls back to 16-mixed / 32). CPU still works via --accelerator cpu.
2. Scale to the full 432,919-nucleus atlas: the expression binary is memory-mapped and read
   per-cell in the DataLoader instead of being fully materialised in RAM. An optional
   --in-memory flag loads everything as float32 when RAM allows (faster).
3. Reproducible log-transform decision (--log-transform auto|on|off) computed from a sampled
   chunk rather than the global max of an in-RAM matrix.
4. CosineAnnealingLR T_max now tracks --epochs (the original hardcoded 40).
5. Everything is argparse-configurable; defaults reproduce the manuscript methodology but on
   GPU and the full atlas.

Input contract (unchanged, produced by the atlas-extraction step — see Code/export_atlas_full.R):
  <atlas_dir>/expr_sub.bin   : flat binary, genes x cells, Fortran order, dtype per --input-dtype
  <atlas_dir>/meta_sub.csv   : one row per cell; columns orig.cluster, orig.ident, dataset
  <atlas_dir>/dims.txt       : "<n_genes> <n_cells>"
  <sig_matrix>               : genes x clusters signature CSV (index_col=0)
"""
import os, sys, json, time, argparse
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping


# --------------------------------------------------------------------------- #
# Model (identical architecture to the CPU draft)
# --------------------------------------------------------------------------- #
class StimulusAutoDecoder(pl.LightningModule):
    """Conditional VAE / auto-decoder for stimulus latent codes."""
    def __init__(self, n_genes, n_clusters, n_organs, n_stages,
                 latent_dim=32, hidden_dim=512, lr=1e-3, epochs=40):
        super().__init__()
        self.save_hyperparameters()
        self.latent_dim = latent_dim
        self.n_genes = n_genes
        self.n_clusters = n_clusters
        self.n_organs = n_organs
        self.n_stages = n_stages

        self.cluster_emb = nn.Embedding(n_clusters, 32)
        self.organ_emb = nn.Embedding(n_organs, 16)
        self.stage_emb = nn.Embedding(n_stages, 16)
        cond_dim = 32 + 16 + 16  # 64

        self.encoder = nn.Sequential(
            nn.Linear(n_genes, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(0.1),
        )
        self.fc_mu = nn.Linear(hidden_dim // 2, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim // 2, latent_dim)

        self.decoder = nn.Sequential(
            nn.Linear(latent_dim + cond_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, n_genes),
        )

        self.classifier = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.GELU(),
            nn.Linear(hidden_dim // 2, n_clusters),
        )

    def encode(self, x):
        h = self.encoder(x)
        return self.fc_mu(h), self.fc_logvar(h)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z, cluster_idx, organ_idx, stage_idx):
        cond = torch.cat([self.cluster_emb(cluster_idx),
                          self.organ_emb(organ_idx),
                          self.stage_emb(stage_idx)], dim=-1)
        return self.decoder(torch.cat([z, cond], dim=-1))

    def forward(self, x, cluster_idx, organ_idx, stage_idx):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z, cluster_idx, organ_idx, stage_idx)
        return recon, mu, logvar, z, self.classifier(z)

    def loss_function(self, recon, x, mu, logvar, class_logits, cluster_idx, beta=0.5):
        recon_loss = F.mse_loss(recon, x, reduction='mean')
        kld = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
        cls_loss = F.cross_entropy(class_logits, cluster_idx)
        return recon_loss + beta * kld + 0.1 * cls_loss, recon_loss, kld, cls_loss

    def training_step(self, batch, batch_idx):
        x, cidx, oidx, sidx = batch
        recon, mu, logvar, z, logits = self(x, cidx, oidx, sidx)
        loss, rl, kld, cls = self.loss_function(recon, x, mu, logvar, logits, cidx)
        self.log_dict({'train_loss': loss, 'recon_loss': rl, 'kld': kld, 'cls_loss': cls},
                      prog_bar=True, on_epoch=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, cidx, oidx, sidx = batch
        recon, mu, logvar, z, logits = self(x, cidx, oidx, sidx)
        loss, rl, kld, cls = self.loss_function(recon, x, mu, logvar, logits, cidx)
        self.log_dict({'val_loss': loss, 'val_recon': rl, 'val_kld': kld, 'val_cls': cls},
                      prog_bar=True, on_epoch=True)
        return loss

    def configure_optimizers(self):
        opt = torch.optim.AdamW(self.parameters(), lr=self.hparams.lr, weight_decay=1e-5)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=self.hparams.epochs)
        return {"optimizer": opt, "lr_scheduler": sched}


# --------------------------------------------------------------------------- #
# Dataset — streams cells from a memory-mapped (genes x cells, Fortran) binary,
# or from an in-RAM float32 array when --in-memory is set.
# --------------------------------------------------------------------------- #
class AtlasDataset(Dataset):
    def __init__(self, source, indices, cluster_idx, organ_idx, stage_idx,
                 log_transform=False, clip_max=10.0, in_memory=False):
        # source: np.memmap (genes x cells) OR np.ndarray (cells x genes, float32)
        self.source = source
        self.indices = indices
        self.in_memory = in_memory
        self.log_transform = log_transform
        self.clip_max = clip_max
        self.cluster_idx = torch.from_numpy(cluster_idx[indices]).long()
        self.organ_idx = torch.from_numpy(organ_idx[indices]).long()
        self.stage_idx = torch.from_numpy(stage_idx[indices]).long()

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, i):
        cell = self.indices[i]
        if self.in_memory:
            vec = self.source[cell]                       # already cells x genes float32
        else:
            vec = np.asarray(self.source[:, cell], dtype=np.float32)  # column = one cell
        if self.log_transform:
            vec = np.log1p(vec)
        np.clip(vec, 0.0, self.clip_max, out=vec)
        return (torch.from_numpy(np.ascontiguousarray(vec, dtype=np.float32)),
                self.cluster_idx[i], self.organ_idx[i], self.stage_idx[i])


# --------------------------------------------------------------------------- #
def decide_log_transform(mode, source, n_cells, in_memory, rng):
    if mode == "on":
        return True
    if mode == "off":
        return False
    # auto: sample up to 20k cells and inspect the max
    k = min(n_cells, 20000)
    sample = rng.choice(n_cells, size=k, replace=False)
    if in_memory:
        mx = float(source[sample].max())
    else:
        mx = float(np.asarray(source[:, sample]).max())
    print(f"[log-transform auto] sampled max = {mx:.3f} -> "
          f"{'log1p' if mx > 20 else 'no transform'}")
    return mx > 20


def main():
    p = argparse.ArgumentParser(description="Full-atlas GPU training of the stimulus auto-decoder.")
    p.add_argument("--atlas-dir", default="/mnt/shared-workspace/autodecoder",
                   help="dir holding expr_sub.bin, meta_sub.csv, dims.txt")
    p.add_argument("--sig-matrix", default="/mnt/shared-workspace/processed/cell_type_signatures.csv")
    p.add_argument("--out-dir", default="/mnt/shared-workspace/autodecoder")
    p.add_argument("--ckpt-dir", default="/workspace/autodecoder_ckpts")
    p.add_argument("--input-dtype", default="float64", choices=["float64", "float32", "float16"],
                   help="dtype of expr_sub.bin (original R writeBin emits float64)")
    p.add_argument("--epochs", type=int, default=40)
    p.add_argument("--batch-size", type=int, default=512)
    p.add_argument("--num-workers", type=int, default=8)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--latent-dim", type=int, default=32)
    p.add_argument("--hidden-dim", type=int, default=512)
    p.add_argument("--accelerator", default="auto", choices=["auto", "gpu", "mps", "cpu"])
    p.add_argument("--devices", default="auto")
    p.add_argument("--precision", default="bf16-mixed",
                   choices=["bf16-mixed", "16-mixed", "32-true", "32"])
    p.add_argument("--log-transform", default="auto", choices=["auto", "on", "off"])
    p.add_argument("--in-memory", action="store_true",
                   help="load the whole atlas into RAM as float32 (faster; needs ~n_cells*n_genes*4 bytes)")
    p.add_argument("--early-stop-patience", type=int, default=10)
    p.add_argument("--seed", type=int, default=42)
    # cardinality guard: expected number of distinct classes (0 = skip that check)
    p.add_argument("--expect-clusters", type=int, default=183)
    p.add_argument("--expect-organs", type=int, default=12)
    p.add_argument("--expect-stages", type=int, default=10)
    p.add_argument("--allow-cardinality-mismatch", action="store_true",
                   help="downgrade a cardinality mismatch from a hard error to a warning")
    args = p.parse_args()

    pl.seed_everything(args.seed, workers=True)
    rng = np.random.default_rng(args.seed)
    torch.set_float32_matmul_precision("high")
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.ckpt_dir, exist_ok=True)

    # --- device resolution + precision guard (CUDA, Apple MPS, or CPU) ---
    mps_ok = getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available()
    accelerator = args.accelerator
    if accelerator == "auto":
        accelerator = "gpu" if torch.cuda.is_available() else ("mps" if mps_ok else "cpu")
    precision = args.precision
    if accelerator == "gpu":
        if precision == "bf16-mixed" and not torch.cuda.is_bf16_supported():
            precision = "16-mixed"
            print("[precision] bf16 unsupported on this GPU -> falling back to 16-mixed")
        print(f"[device] CUDA: {torch.cuda.get_device_name(0)} "
              f"({torch.cuda.get_device_properties(0).total_memory/1e9:.1f} GB)")
    elif accelerator == "mps":
        # Apple Metal: bf16 autocast isn't supported; default to full precision.
        if precision == "bf16-mixed":
            precision = "32-true"
            print("[precision] Apple MPS -> 32-true (bf16 not supported on Metal)")
        print("[device] Apple GPU (MPS / Metal). Tip: use --in-memory (avoids memmap in "
              "DataLoader workers); --num-workers 0 is safest on macOS.")
    else:  # cpu
        precision = "32-true"
    print(f"[device] accelerator={accelerator} devices={args.devices} precision={precision}")

    # --- reference / signature matrix ---
    print("=== Loading atlas reference data ===")
    sig_df = pd.read_csv(args.sig_matrix, index_col=0)
    hvgs = list(sig_df.index)
    print(f"Signature matrix: {sig_df.shape} (genes x clusters)")

    # --- dims + metadata ---
    with open(os.path.join(args.atlas_dir, "dims.txt")) as f:
        n_genes_file, n_cells_file = (int(x) for x in f.read().split())
    meta_sub = pd.read_csv(os.path.join(args.atlas_dir, "meta_sub.csv"))
    n_cells = len(meta_sub)
    assert n_cells == n_cells_file, f"Cell count mismatch: {n_cells} vs {n_cells_file}"
    n_hvgs = n_genes_file
    print(f"Cells: {n_cells:,}  Genes: {n_hvgs:,}  "
          f"(expecting the FULL atlas — first draft used ~60k)")

    # --- labels ---
    # Required-column guard: fail loudly (not a cryptic KeyError) if the export schema is wrong.
    required = {"cluster": "orig.cluster", "organ": "orig.ident", "stage": "dataset"}
    missing_cols = [c for c in required.values() if c not in meta_sub.columns]
    if missing_cols:
        raise SystemExit(
            f"meta_sub.csv is missing required column(s): {missing_cols}. "
            f"Columns present: {list(meta_sub.columns)}. "
            f"Re-export with Code/export_atlas.R (it emits orig.cluster / orig.ident / dataset).")

    cluster_labels = meta_sub['orig.cluster'].astype(str).values
    organ_labels = meta_sub['orig.ident'].astype(str).values
    stage_labels = meta_sub['dataset'].astype(str).values
    cluster_to_idx = {c: i for i, c in enumerate(sorted(set(cluster_labels)))}
    organ_to_idx = {o: i for i, o in enumerate(sorted(set(organ_labels)))}
    stage_to_idx = {s: i for i, s in enumerate(sorted(set(stage_labels)))}
    cluster_idx = np.array([cluster_to_idx[c] for c in cluster_labels])
    organ_idx = np.array([organ_to_idx[o] for o in organ_labels])
    stage_idx = np.array([stage_to_idx[s] for s in stage_labels])
    n_cluster_classes, n_organ_classes, n_stage_classes = \
        len(cluster_to_idx), len(organ_to_idx), len(stage_to_idx)
    print(f"Cluster classes: {n_cluster_classes}, Organ classes: {n_organ_classes}, "
          f"Stage classes: {n_stage_classes}")

    # Cardinality guard: the exported metadata must match the atlas structure the
    # architecture assumes (embeddings sized for 183 clusters / 12 organs / 10 stages).
    # A mismatch usually means a wrong object, wrong column mapping, or a partial export —
    # block training by default so it can't slip through unnoticed.
    card_problems = []
    for name, got, exp in (("clusters", n_cluster_classes, args.expect_clusters),
                           ("organs", n_organ_classes, args.expect_organs),
                           ("stages", n_stage_classes, args.expect_stages)):
        if exp and got != exp:
            card_problems.append(f"{name}: found {got}, expected {exp}")
    if card_problems:
        msg = ("Cardinality guard: " + "; ".join(card_problems) +
               ". Check the atlas object / column mappings in Code/export_atlas.R, "
               "or set --expect-* to the correct values (0 disables a check).")
        if args.allow_cardinality_mismatch:
            print("WARNING: " + msg)
        else:
            raise SystemExit("ERROR: " + msg + " Re-export, or pass --allow-cardinality-mismatch "
                             "to proceed anyway (embeddings will size to the data).")
    else:
        print("[cardinality guard] OK — clusters/organs/stages match expectations.")

    # --- expression source (memmap, or in-RAM float32) ---
    dtype = np.dtype(args.input_dtype)
    bin_path = os.path.join(args.atlas_dir, "expr_sub.bin")
    expr_mm = np.memmap(bin_path, dtype=dtype, mode="r",
                        shape=(n_hvgs, n_cells), order="F")  # genes x cells
    gb = n_hvgs * n_cells * 4 / 1e9
    if args.in_memory:
        print(f"[data] loading full atlas into RAM as float32 (~{gb:.1f} GB) ...")
        source = np.ascontiguousarray(expr_mm.T, dtype=np.float32)  # cells x genes
        in_memory = True
    else:
        print(f"[data] streaming cells from memmap (would be ~{gb:.1f} GB in RAM if loaded)")
        source = expr_mm
        in_memory = False

    do_log = decide_log_transform(args.log_transform, source, n_cells, in_memory, rng)

    # --- train/val split ---
    perm = rng.permutation(n_cells)
    n_val = max(1000, int(0.1 * n_cells))
    val_idx, train_idx = perm[:n_val], perm[n_val:]
    common = dict(log_transform=do_log, in_memory=in_memory)
    train_ds = AtlasDataset(source, train_idx, cluster_idx, organ_idx, stage_idx, **common)
    val_ds = AtlasDataset(source, val_idx, cluster_idx, organ_idx, stage_idx, **common)
    print(f"Train: {len(train_ds):,}  Val: {len(val_ds):,}")

    pin = accelerator == "gpu"
    loader_kw = dict(batch_size=args.batch_size, num_workers=args.num_workers,
                     pin_memory=pin, persistent_workers=args.num_workers > 0)
    train_loader = DataLoader(train_ds, shuffle=True, **loader_kw)
    val_loader = DataLoader(val_ds, shuffle=False, **loader_kw)

    # --- model ---
    model = StimulusAutoDecoder(
        n_genes=n_hvgs, n_clusters=n_cluster_classes, n_organs=n_organ_classes,
        n_stages=n_stage_classes, latent_dim=args.latent_dim,
        hidden_dim=args.hidden_dim, lr=args.lr, epochs=args.epochs)
    print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")

    checkpoint_cb = ModelCheckpoint(
        dirpath=args.ckpt_dir, filename="autodecoder-{epoch:02d}-{val_loss:.4f}",
        save_top_k=3, monitor="val_loss", mode="min")
    early_stop_cb = EarlyStopping(monitor="val_loss", patience=args.early_stop_patience, mode="min")

    trainer = pl.Trainer(
        max_epochs=args.epochs,
        accelerator=accelerator,
        devices=args.devices,
        precision=precision,
        callbacks=[checkpoint_cb, early_stop_cb],
        log_every_n_steps=20,
        enable_progress_bar=True,
        default_root_dir=os.path.join(args.ckpt_dir, "lightning"),
        logger=pl.loggers.CSVLogger(save_dir=args.ckpt_dir, name="autodecoder"))

    print("\n=== Training auto-decoder (full atlas, GPU) ===")
    t0 = time.time()
    trainer.fit(model, train_loader, val_loader)
    print(f"\nTraining completed in {(time.time()-t0)/60:.1f} min")
    best_path = checkpoint_cb.best_model_path
    print(f"Best model: {best_path}")
    if checkpoint_cb.best_model_score is not None:
        print(f"Best val_loss: {float(checkpoint_cb.best_model_score):.4f}")

    with open(os.path.join(args.out_dir, "model_encodings.json"), "w") as f:
        json.dump({'cluster_to_idx': cluster_to_idx, 'organ_to_idx': organ_to_idx,
                   'stage_to_idx': stage_to_idx, 'hvgs': hvgs, 'latent_dim': args.latent_dim,
                   'n_hvgs': n_hvgs, 'best_model_path': best_path,
                   'n_cells_trained': int(n_cells)}, f, indent=2)

    # --- atlas-wide latent embeddings (batched on device) ---
    print("\n=== Computing atlas latent embeddings ===")
    device = torch.device("cuda" if (accelerator == "gpu" and torch.cuda.is_available())
                          else "mps" if accelerator == "mps" else "cpu")
    model.to(device).eval()
    all_ds = AtlasDataset(source, np.arange(n_cells), cluster_idx, organ_idx, stage_idx, **common)
    all_loader = DataLoader(all_ds, batch_size=max(1024, args.batch_size),
                            shuffle=False, num_workers=args.num_workers, pin_memory=pin)
    zs, preds = [], []
    with torch.no_grad():
        for x, cidx, oidx, sidx in all_loader:
            x = x.to(device, non_blocking=True)
            mu, _ = model.encode(x)
            zs.append(mu.float().cpu().numpy())
            preds.append(model.classifier(mu).argmax(dim=-1).cpu().numpy())
    all_z = np.vstack(zs)
    all_cls_pred = np.concatenate(preds)
    print(f"Latent embeddings: {all_z.shape}")

    np.save(os.path.join(args.out_dir, "atlas_latent_embeddings.npy"), all_z)
    np.save(os.path.join(args.out_dir, "atlas_cls_predictions.npy"), all_cls_pred)

    # per-cluster mean latent = stimulus code per cell type (also the web tool's Phase-2 artifact)
    idx_to_cluster = {v: k for k, v in cluster_to_idx.items()}
    cluster_latent = {}
    for cl in range(n_cluster_classes):
        mask = cluster_idx == cl
        if mask.sum() > 0:
            cluster_latent[idx_to_cluster[cl]] = all_z[mask].mean(axis=0).tolist()
    with open(os.path.join(args.out_dir, "cluster_stimulus_codes.json"), "w") as f:
        json.dump(cluster_latent, f, indent=2)

    print("\n=== DONE Subsystem 3 (GPU training) ===")
    print(f"Model:              {best_path}")
    print(f"Latent embeddings:  {args.out_dir}/atlas_latent_embeddings.npy")
    print(f"Stimulus codes:     {args.out_dir}/cluster_stimulus_codes.json")


if __name__ == "__main__":
    main()
