'''
Model part for LCL 
2025-4-11
'''

# general package
import tempfile
import os
import numpy as np

# deep learning package
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.multiprocessing import cpu_count
import torchvision.transforms as T


#-----------------------------------------Add projection Head for embedding----------------------------------------------
###NEW

# Transformer-based encoder tailored for single-cell contrastive learning
class TransformerCellEncoder(nn.Module):
    def __init__(self, num_genes=2000, embed_dim=128, num_heads=8, hidden_dim=256, final_dim=64):
        super(TransformerCellEncoder, self).__init__()
        self.input_embedding = nn.Linear(num_genes, embed_dim)
        self.mha = nn.MultiheadAttention(embed_dim, num_heads, batch_first=True)
        self.norm1 = nn.LayerNorm(embed_dim)
        self.ff = nn.Sequential(
            nn.Linear(embed_dim, hidden_dim),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_dim, embed_dim),
        )
        self.norm2 = nn.LayerNorm(embed_dim)
        self.final_linear = nn.Linear(embed_dim, final_dim)

    def forward(self, x):
        """
        Input shape: (batch_size, cells_per_lineage, num_genes)
        Output shape: (batch_size, cells_per_lineage, final_dim)
        """
        x_emb = self.input_embedding(x)                             # (B, C, embed_dim)
        assert x_emb.dim() == 3, f"x_emb must be 3D (B, C, embed_dim), but got shape {x_emb.shape}"
        attn_out, _ = self.mha(x_emb, x_emb, x_emb)                 # (B, C, embed_dim)
        x = self.norm1(x_emb + attn_out)                            # (B, C, embed_dim)
        ff_out = self.ff(x)                                         # (B, C, embed_dim)
        x = self.norm2(x + ff_out)                                  # (B, C, embed_dim)
        final_embedding = self.final_linear(x)                      # (B, C, final_dim)
        return final_embedding

# Projection head integrated with transformer encoder for contrastive learning
class AddProjectionTransformer(nn.Module):
    def __init__(self, config):
        """
        config attributes:
            - input_dim (int): Number of input genes (features).
            - embed_dim (int): Embedding dimension after initial projection.
            - num_heads (int): Number of attention heads.
            - hidden_dim (int): Hidden dimension in feed-forward layers.
            - final_dim (int): Output dimension from Transformer encoder.
            - projection_hidden_dim (int, optional): Hidden dimension in projection MLP.
            - embedding_size (int): Final output embedding dimension.
        """
        super(AddProjectionTransformer, self).__init__()
        self.base_encoder = TransformerCellEncoder(
            num_genes=config.num_genes,
            embed_dim=getattr(config, 'embed_dim', 128),
            num_heads=getattr(config, 'num_heads', 8),
            hidden_dim=getattr(config, 'hidden_dim', 256),
            final_dim=getattr(config, 'final_dim', 64)
        )

        proj_hidden_dim = getattr(config, 'projection_hidden_dim', config.final_dim)
        self.projection = nn.Sequential(
            nn.Linear(config.final_dim, proj_hidden_dim),
            nn.BatchNorm1d(proj_hidden_dim),
            nn.ReLU(inplace=True),
            nn.Linear(proj_hidden_dim, config.embedding_size),
            nn.BatchNorm1d(config.embedding_size),
        )

    def forward(self, x, return_embedding=False):
        """
        x shape: (batch_size, cells_per_lineage, input_dim)
        - If return_embedding=True, returns base encoder features: (batch_size * cells_per_lineage, final_dim)
        - Otherwise, returns projected embeddings: (batch_size * cells_per_lineage, embedding_size)
        """
        base_embedding = self.base_encoder(x)  # (B, C, final_dim)

        # Flatten batch and cells_per_lineage dimensions
        B, C, final_dim = base_embedding.shape
        flat_embedding = base_embedding.view(B * C, final_dim)

        if return_embedding:
            return flat_embedding  # (B * C, final_dim)

        proj_embedding = self.projection(flat_embedding)  # (B * C, embedding_size)

        # print("Embedding mean:", proj_embedding.mean().item())
        # print("Embedding std:", proj_embedding.std().item())


        return proj_embedding

    def get_features(self, x):
        """
        Extract base encoder features directly.
        x shape: (batch_size, cells_per_lineage, input_dim)
        output shape: (batch_size * cells_per_lineage, final_dim)
        """
        with torch.no_grad():
            base_embedding = self.base_encoder(x)
            B, C, final_dim = base_embedding.shape
            return base_embedding.view(B * C, final_dim)




#--------------------------------------------------------ContrastiveLoss-------------------------------------------------
###NEW


class ContrastiveLoss(nn.Module):
    def __init__(self, temperature=0.5):
        """
        Contrastive Loss for supervised contrastive learning on cell embeddings.

        Args:
            temperature (float): Scaling factor for similarity scores.
        """
        super(ContrastiveLoss, self).__init__()
        self.temperature = temperature

    def forward(self, features):
        """
        Compute the supervised contrastive loss for the provided features.

        Args:
            features (torch.Tensor): Embeddings of shape
                                     (batch_size, cells_per_lineage, embedding_dim).

        Returns:
            torch.Tensor: Scalar contrastive loss.
        """
        batch_size, cells_per_lineage, embedding_dim = features.shape

        # NaN check before anything
        if torch.isnan(features).any():
            print("[DEBUG] NaN detected in input features!")
            print("features:", features)
            exit()

        # Normalize features
        features = F.normalize(features, dim=2)

        # Reshape to (batch_size * cells_per_lineage, embedding_dim)
        features_flat = features.view(batch_size * cells_per_lineage, embedding_dim)

        # Compute similarity matrix
        similarity_matrix = torch.matmul(features_flat, features_flat.T) / self.temperature

        # Create labels indicating positive pairs (cells from same lineage)
        labels = torch.arange(batch_size).repeat_interleave(cells_per_lineage).to(features.device)

        # Mask to exclude self-comparisons
        mask = torch.eye(batch_size * cells_per_lineage, dtype=torch.bool).to(features.device)

        # Compute log-softmax of similarities
        logits = similarity_matrix.masked_fill(mask, float('-inf'))
        log_prob = F.log_softmax(logits, dim=1)

        # Create mask for positive pairs (same lineage but not the same cell)
        positive_mask = labels.unsqueeze(0) == labels.unsqueeze(1)
        positive_mask = positive_mask & (~mask)  # remove self-pairs

        # Count positive pairs for normalization
        positive_count = positive_mask.sum(1)

        # Compute loss (only for cells with at least one positive pair)
        masked_log_prob = log_prob.masked_fill(~positive_mask, 0.0)
        loss = -masked_log_prob.sum(1) / positive_count.clamp(min=1)

        loss = loss.mean()

        # Post-computation check
        if torch.isnan(loss) or torch.isinf(loss):
            print("[DEBUG] NaN or Inf detected in loss!")
            print("loss:", loss)
            print("features_flat (sample):", features_flat[0])
            print("similarity_matrix (sample):", similarity_matrix[0][:10])
            print("log_prob (sample):", log_prob[0][:10])
            exit()

        return loss