#!/usr/bin/env python3
"""Subsystem 3: Custom stimulus auto-decoder (conditional VAE).
Trains a variational auto-decoder on the Salk atlas to learn cell-type- and
stimulus-specific latent codes, then projects bulk samples onto the learned
latent space to recover cell-type-specific tropism response scores.

Architecture:
- Encoder: cell expression vector (HVGs) -> latent z (stimulus code, dim 32)
- Decoder: (z, cell_type_code, dev_stage_code) -> reconstructed expression
- Auxiliary head: cell-type classification from z
- KL regularization on z toward Gaussian prior

Since no GPU is available, we subsample the atlas to ~60k nuclei covering all
major cell types and train for ~40 epochs on CPU. The full cell-type signature
matrix (183 clusters x 4000 HVGs) is still used for deconvolution.
"""
import os, sys, json, time, random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping

# Paths
ATLAS_COMPACT = "/mnt/shared-workspace/processed/atlas_compact_reference.rds"
SIG_MATRIX = "/mnt/shared-workspace/processed/cell_type_signatures.csv"
CLUSTER_SUMMARY = "/mnt/shared-workspace/processed/cluster_summary.csv"
OUT_DIR = "/mnt/shared-workspace/autodecoder"
os.makedirs(OUT_DIR, exist_ok=True)

# Load the signature matrix and cluster info
print("=== Loading atlas reference data ===")
sig_df = pd.read_csv(SIG_MATRIX, index_col=0)
print(f"Signature matrix: {sig_df.shape} (genes x clusters)")
clusters = list(sig_df.columns)
n_clusters = len(clusters)
n_hvgs = sig_df.shape[0]
hvgs = list(sig_df.index)

cluster_summary = pd.read_csv(CLUSTER_SUMMARY)
print(f"Cluster summary: {cluster_summary.shape}")
print(f"Clusters: {n_clusters}, HVGs: {n_hvgs}")

# Load the pre-extracted per-cell HVG expression (from fix_signatures.R)
print("\n=== Loading pre-extracted expression ===")
meta_sub = pd.read_csv("/mnt/shared-workspace/autodecoder/meta_sub.csv")
n_cells = len(meta_sub)
# Read dimensions
with open("/mnt/shared-workspace/autodecoder/dims.txt") as f:
    dims = f.read().split()
n_genes_file, n_cells_file = int(dims[0]), int(dims[1])
assert n_cells == n_cells_file, f"Cell count mismatch: {n_cells} vs {n_cells_file}"
n_hvgs = n_genes_file
print(f"Cells: {n_cells}, Genes: {n_hvgs}")

# Read the binary expression matrix (column-major / Fortran order, as saved by R writeBin)
expr_flat = np.fromfile("/mnt/shared-workspace/autodecoder/expr_sub.bin", dtype=np.float64)
expr = expr_flat.reshape(n_hvgs, n_cells, order='F')  # genes x cells
print(f"Expression matrix: {expr.shape} (genes x cells)")
print(f"Expression range: [{expr.min():.3f}, {expr.max():.3f}]")
print(f"Mean: {expr.mean():.3f}")
del expr_flat

# Prepare labels
cluster_labels = meta_sub['orig.cluster'].astype(str).values
organ_labels = meta_sub['orig.ident'].astype(str).values
stage_labels = meta_sub['dataset'].astype(str).values

# Encode labels
cluster_to_idx = {c: i for i, c in enumerate(sorted(set(cluster_labels)))}
organ_to_idx = {o: i for i, o in enumerate(sorted(set(organ_labels)))}
stage_to_idx = {s: i for i, s in enumerate(sorted(set(stage_labels)))}

cluster_idx = np.array([cluster_to_idx[c] for c in cluster_labels])
organ_idx = np.array([organ_to_idx[o] for o in organ_labels])
stage_idx = np.array([stage_to_idx[s] for s in stage_labels])

n_cluster_classes = len(cluster_to_idx)
n_organ_classes = len(organ_to_idx)
n_stage_classes = len(stage_to_idx)
print(f"Cluster classes: {n_cluster_classes}, Organ classes: {n_organ_classes}, Stage classes: {n_stage_classes}")

# Transpose expression for PyTorch (cells x genes)
expr_T = expr.T.astype(np.float32)  # cells x genes
del expr

# Log-transform if not already (Seurat data layer is log-normalized)
# Check if values look log-scale (max < 15) or count-scale
if expr_T.max() > 20:
    print("Applying log1p transform...")
    expr_T = np.log1p(expr_T)

# Clip extreme values
expr_T = np.clip(expr_T, 0, 10)

print(f"Final training matrix: {expr_T.shape}, range [{expr_T.min():.3f}, {expr_T.max():.3f}]")

# === Model Definition ===
class StimulusAutoDecoder(pl.LightningModule):
    """Conditional VAE / auto-decoder for stimulus latent codes."""
    def __init__(self, n_genes, n_clusters, n_organs, n_stages,
                 latent_dim=32, hidden_dim=512, lr=1e-3):
        super().__init__()
        self.save_hyperparameters()
        self.latent_dim = latent_dim
        self.n_genes = n_genes
        self.n_clusters = n_clusters
        self.n_organs = n_organs
        self.n_stages = n_stages

        # Embeddings for conditioning
        self.cluster_emb = nn.Embedding(n_clusters, 32)
        self.organ_emb = nn.Embedding(n_organs, 16)
        self.stage_emb = nn.Embedding(n_stages, 16)
        cond_dim = 32 + 16 + 16  # 64

        # Encoder: expression -> z params
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

        # Decoder: (z, conditioning) -> expression
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

        # Auxiliary cell-type classification head
        self.classifier = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.GELU(),
            nn.Linear(hidden_dim // 2, n_clusters),
        )

    def encode(self, x):
        h = self.encoder(x)
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        return mu, logvar

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z, cluster_idx, organ_idx, stage_idx):
        c_emb = self.cluster_emb(cluster_idx)
        o_emb = self.organ_emb(organ_idx)
        s_emb = self.stage_emb(stage_idx)
        cond = torch.cat([c_emb, o_emb, s_emb], dim=-1)
        z_cond = torch.cat([z, cond], dim=-1)
        return self.decoder(z_cond)

    def forward(self, x, cluster_idx, organ_idx, stage_idx):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z, cluster_idx, organ_idx, stage_idx)
        class_logits = self.classifier(z)
        return recon, mu, logvar, z, class_logits

    def loss_function(self, recon, x, mu, logvar, class_logits, cluster_idx, beta=0.5):
        # Reconstruction loss (MSE)
        recon_loss = F.mse_loss(recon, x, reduction='mean')
        # KL divergence
        kld = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
        # Classification loss (auxiliary)
        cls_loss = F.cross_entropy(class_logits, cluster_idx)
        return recon_loss + beta * kld + 0.1 * cls_loss, recon_loss, kld, cls_loss

    def training_step(self, batch, batch_idx):
        x, cidx, oidx, sidx = batch
        recon, mu, logvar, z, class_logits = self(x, cidx, oidx, sidx)
        loss, rl, kld, cls = self.loss_function(recon, x, mu, logvar, class_logits, cidx)
        self.log_dict({'train_loss': loss, 'recon_loss': rl, 'kld': kld, 'cls_loss': cls},
                      prog_bar=True, on_epoch=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, cidx, oidx, sidx = batch
        recon, mu, logvar, z, class_logits = self(x, cidx, oidx, sidx)
        loss, rl, kld, cls = self.loss_function(recon, x, mu, logvar, class_logits, cidx)
        self.log_dict({'val_loss': loss, 'val_recon': rl, 'val_kld': kld, 'val_cls': cls},
                      prog_bar=True, on_epoch=True)
        return loss

    def configure_optimizers(self):
        opt = torch.optim.AdamW(self.parameters(), lr=self.hparams.lr, weight_decay=1e-5)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=40)
        return {"optimizer": opt, "lr_scheduler": sched}


class AtlasDataset(Dataset):
    def __init__(self, expr, cluster_idx, organ_idx, stage_idx):
        self.expr = torch.from_numpy(expr)
        self.cluster_idx = torch.from_numpy(cluster_idx).long()
        self.organ_idx = torch.from_numpy(organ_idx).long()
        self.stage_idx = torch.from_numpy(stage_idx).long()

    def __len__(self):
        return self.expr.shape[0]

    def __getitem__(self, idx):
        return self.expr[idx], self.cluster_idx[idx], self.organ_idx[idx], self.stage_idx[idx]


# === Training ===
print("\n=== Setting up training ===")
# Split into train/val
np.random.seed(42)
n_total = expr_T.shape[0]
perm = np.random.permutation(n_total)
n_val = max(1000, int(0.1 * n_total))
val_idx = perm[:n_val]
train_idx = perm[n_val:]

train_ds = AtlasDataset(expr_T[train_idx], cluster_idx[train_idx], organ_idx[train_idx], stage_idx[train_idx])
val_ds = AtlasDataset(expr_T[val_idx], cluster_idx[val_idx], organ_idx[val_idx], stage_idx[val_idx])

print(f"Train: {len(train_ds)}, Val: {len(val_ds)}")

train_loader = DataLoader(train_ds, batch_size=256, shuffle=True, num_workers=2, persistent_workers=True)
val_loader = DataLoader(val_ds, batch_size=256, shuffle=False, num_workers=2, persistent_workers=True)

# Model
model = StimulusAutoDecoder(
    n_genes=n_hvgs,
    n_clusters=n_cluster_classes,
    n_organs=n_organ_classes,
    n_stages=n_stage_classes,
    latent_dim=32,
    hidden_dim=512,
    lr=1e-3
)
print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")

# Callbacks
os.makedirs("/workspace/autodecoder_ckpts", exist_ok=True)
checkpoint_cb = ModelCheckpoint(
    dirpath="/workspace/autodecoder_ckpts", filename="autodecoder-{epoch:02d}-{val_loss:.4f}",
    save_top_k=3, monitor="val_loss", mode="min"
)
early_stop_cb = EarlyStopping(monitor="val_loss", patience=8, mode="min")

# Trainer
trainer = pl.Trainer(
    max_epochs=40,
    accelerator="cpu",
    devices=1,
    callbacks=[checkpoint_cb, early_stop_cb],
    log_every_n_steps=20,
    enable_progress_bar=True,
    default_root_dir="/workspace/autodecoder_lightning",
    logger=pl.loggers.CSVLogger(save_dir="/workspace/autodecoder_lightning", name="autodecoder"),
)

print("\n=== Training auto-decoder ===")
t0 = time.time()
trainer.fit(model, train_loader, val_loader)
t1 = time.time()
print(f"\nTraining completed in {(t1-t0)/60:.1f} min")

# Save final model
best_path = checkpoint_cb.best_model_path
print(f"Best model: {best_path}")
print(f"Best val_loss: {checkpoint_cb.best_model_score:.4f}")

# Save label encodings for downstream use
encodings = {
    'cluster_to_idx': cluster_to_idx,
    'organ_to_idx': organ_to_idx,
    'stage_to_idx': stage_to_idx,
    'hvgs': hvgs,
    'latent_dim': 32,
    'n_hvgs': n_hvgs,
    'best_model_path': best_path,
}
with open(f"{OUT_DIR}/model_encodings.json", "w") as f:
    json.dump(encodings, f, indent=2)

# === Compute atlas-wide latent embeddings ===
print("\n=== Computing atlas latent embeddings ===")
model.eval()
all_loader = DataLoader(AtlasDataset(expr_T, cluster_idx, organ_idx, stage_idx),
                        batch_size=512, shuffle=False, num_workers=2)
all_z = []
all_cls_pred = []
with torch.no_grad():
    for batch in all_loader:
        x, cidx, oidx, sidx = batch
        mu, logvar = model.encode(x)
        z = mu  # use mean for inference
        class_logits = model.classifier(z)
        all_z.append(mu.numpy())
        all_cls_pred.append(class_logits.argmax(dim=-1).numpy())

all_z = np.vstack(all_z)
all_cls_pred = np.concatenate(all_cls_pred)
print(f"Latent embeddings: {all_z.shape}")

# Save embeddings
np.save(f"{OUT_DIR}/atlas_latent_embeddings.npy", all_z)
np.save(f"{OUT_DIR}/atlas_cls_predictions.npy", all_cls_pred)

# Per-cluster mean latent (stimulus code per cell type)
cluster_latent = {}
for cl_idx in range(n_cluster_classes):
    mask = cluster_idx == cl_idx
    if mask.sum() > 0:
        cl_name = [k for k, v in cluster_to_idx.items() if v == cl_idx][0]
        cluster_latent[cl_name] = all_z[mask].mean(axis=0).tolist()

with open(f"{OUT_DIR}/cluster_stimulus_codes.json", "w") as f:
    json.dump(cluster_latent, f, indent=2)

print(f"\n=== DONE Subsystem 3 (training) ===")
print(f"Model saved to: {best_path}")
print(f"Latent embeddings: {OUT_DIR}/atlas_latent_embeddings.npy")
print(f"Cluster stimulus codes: {OUT_DIR}/cluster_stimulus_codes.json")
