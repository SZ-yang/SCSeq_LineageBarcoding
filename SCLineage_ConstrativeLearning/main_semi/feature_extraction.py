#!/usr/bin/env python
import os
import time
import argparse
import numpy as np
import anndata as ad
import torch
from torch.utils.data import DataLoader, TensorDataset

from LCL_Model_Semi import BaseEncoder_ProjHead_MLP
# ContrastiveLoss is not needed here unless you plan to re-train
# from LCL_Model_Semi import ContrastiveLoss


def get_args():
    parser = argparse.ArgumentParser(
        description="Extract base-encoder features from a trained LCL model"
    )
    parser.add_argument(
        "--inputFilePath", required=True,
        help="Path to AnnData (.h5ad) for feature extraction"
    )
    parser.add_argument(
        "--batch_size", type=int, default=50,
        help="Batch size for featurizing all cells"
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Directory to save extracted feature .npy"
    )
    parser.add_argument(
        "--resume_from_checkpoint", required=True,
        help="Path to trained .ckpt file"
    )
    parser.add_argument(
        "--hidden_dims", default="1024,256,64",
        help="Comma-separated hidden dimensions for encoder"
    )
    parser.add_argument(
        "--embedding_size", type=int, default=32,
        help="Projection head output dimension"
    )
    return parser.parse_args()


class Hparams:
    def __init__(self, args):
        self.input_dim = None  # will be set from AnnData
        self.hidden_dims = [int(x) for x in args.hidden_dims.split(",")]
        self.embedding_size = args.embedding_size

        # these follow your training defaults (unused here)
        self.unlabeled_per_batch = 5
        self.lambda_penalty = 1.0
        self.epochs = 220
        self.batch_size = args.batch_size
        self.size_factor = 0.3
        self.temperature = 0.5

        # I/O
        self.out_dir = args.output_dir
        self.ckpt_path = args.resume_from_checkpoint
        self.file_path = args.inputFilePath


def load_checkpoint(model, ckpt_path):
    checkpoint = torch.load(ckpt_path, map_location="cpu")
    model.load_state_dict(checkpoint["state_dict"], strict=False)
    print(f"Loaded checkpoint from {ckpt_path}")
    return model


def main():
    start_time = time.time()
    args = get_args()
    cfg = Hparams(args)

    # 1) infer real input_dim
    adata = ad.read_h5ad(cfg.file_path)
    cfg.input_dim = adata.n_vars
    print(f"Detected input_dim = {cfg.input_dim}")

    # 2) build model and load weights
    model = BaseEncoder_ProjHead_MLP(cfg)
    model = load_checkpoint(model, cfg.ckpt_path)
    model.eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(device)

    # 3) prepare count matrix → float32 tensor
    X = adata.X
    try:
        X = X.toarray()
    except AttributeError:
        pass
    X = np.asarray(X, dtype=np.float32)
    tX = torch.from_numpy(X)

    ds = TensorDataset(tX)
    loader = DataLoader(
        ds,
        batch_size=cfg.batch_size,
        shuffle=False,
        num_workers=1,
        drop_last=False,
    )
    print(f"Featurizing {adata.n_obs} cells in {len(loader)} batches")

    # 4) extract features
    feats = []
    with torch.no_grad():
        for batch in loader:
            x = batch[0].to(device)
            z = model.get_features(x)
            feats.append(z.cpu().numpy())

    feats = np.concatenate(feats, axis=0)
    print("Extracted features shape:", feats.shape)

    # 5) save
    os.makedirs(cfg.out_dir, exist_ok=True)
    out_path = os.path.join(
        cfg.out_dir,
        f"scBaseEncoderFeat_Z_bs{cfg.batch_size}_tau{cfg.temperature}.npy"
    )
    np.save(out_path, feats)
    print(f"Saved features to {out_path}")

    elapsed = (time.time() - start_time) / 60
    print(f"Done in {elapsed:.2f} minutes")


if __name__ == "__main__":
    main()