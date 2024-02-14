# general package
import tempfile
import os
import scipy
import numpy as np
import pandas as pd
import copy
import random

# single cell package
import scanpy as sc
import anndata as ad

# deep learning package
import torch
import torchvision.models as models
import numpy as np
import os
import torchvision.transforms as T
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.multiprocessing import cpu_count
import torchvision.transforms as T


def default(val, def_val):
    return def_val if val is None else val

def reproducibility(config):
    SEED = int(config.seed)
    torch.manual_seed(SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    np.random.seed(SEED)
    if (config.cuda):
        torch.cuda.manual_seed(SEED)


def device_as(t1, t2):
    """
    Moves t1 to the device of t2
    """
    return t1.to(t2.device)

# From https://github.com/PyTorchLightning/pytorch-lightning/issues/924
def weights_update(model, checkpoint_path):
    checkpoint = torch.load(checkpoint_path)
    model_dict = model.state_dict()
    pretrained_dict = {k: v for k, v in checkpoint['state_dict'].items() if k in model_dict}
    model_dict.update(pretrained_dict)
    model.load_state_dict(model_dict)
    print(f'Checkpoint {checkpoint_path} was loaded')
    return model



class ContrastiveLoss(nn.Module):
    """
    Vanilla Contrastive loss, also called InfoNceLoss as in SimCLR paper
    """
    def __init__(self, batch_size, temperature=0.5):
        super().__init__()
        self.batch_size = batch_size
        self.temperature = temperature
        self.mask = (~torch.eye(batch_size * 2, batch_size * 2, dtype=bool)).float()

    def calc_similarity_batch(self, a, b):
        representations = torch.cat([a, b], dim=0)
        return F.cosine_similarity(representations.unsqueeze(1), representations.unsqueeze(0), dim=2)

    def forward(self, proj_1, proj_2):
        """
        proj_1 and proj_2 are batched embeddings [batch, embedding_dim]
        where corresponding indices are pairs
        z_i, z_j in the SimCLR paper
        """
        batch_size = proj_1.shape[0]
        z_i = F.normalize(proj_1, p=2, dim=1)
        z_j = F.normalize(proj_2, p=2, dim=1)

        similarity_matrix = self.calc_similarity_batch(z_i, z_j)

        sim_ij = torch.diag(similarity_matrix, batch_size)
        sim_ji = torch.diag(similarity_matrix, -batch_size)

        positives = torch.cat([sim_ij, sim_ji], dim=0)

        nominator = torch.exp(positives / self.temperature)

        denominator = device_as(self.mask, similarity_matrix) * torch.exp(similarity_matrix / self.temperature)

        all_losses = -torch.log(nominator / torch.sum(denominator, dim=1))
        loss = torch.sum(all_losses) / (2 * self.batch_size)
        return loss


#-------------------------------------------------------------------------------------------------------------------
"""## Add projection Head for embedding and training logic with pytorch lightning model"""
### 改掉resnet  换成mlp

import pytorch_lightning as pl
import torch.nn.functional as F
from pl_bolts.optimizers.lr_scheduler import LinearWarmupCosineAnnealingLR
from torch.optim import SGD, Adam


class AddProjectionMLP(nn.Module):
    

    def __init__(self, config):

        """
        input_dim: The size of the input features. 
        hidden_dims: A list of integers where each integer specifies the number of neurons in that hidden layer.
        embedding_size: The size of the output from the projection head, which is also the size of the final embeddings.
        """
        super(AddProjectionMLP, self).__init__()

        # Define the MLP as the base encoder
        input_dim = config.input_dim
        hidden_dims = config.hidden_dims
        embedding_size = config.embedding_size

        layers = []
        for i in range(len(hidden_dims)):
            layers.append(nn.Linear(input_dim if i == 0 else hidden_dims[i-1], hidden_dims[i]))
            layers.append(nn.BatchNorm1d(hidden_dims[i]))
            layers.append(nn.ReLU(inplace=True))
        self.base_encoder = nn.Sequential(*layers)

        # Define the projection head
        self.projection = nn.Sequential(
            nn.Linear(in_features=hidden_dims[-1], out_features=hidden_dims[-1]),
            nn.BatchNorm1d(hidden_dims[-1]),
            nn.ReLU(),
            nn.Linear(in_features=hidden_dims[-1], out_features=embedding_size),
            nn.BatchNorm1d(embedding_size),
        )

    def forward(self, x, return_embedding=False):
        # Flatten the input if necessary
        x = x.view(x.size(0), -1)
        embedding = self.base_encoder(x)
        if return_embedding:
            return embedding
        return self.projection(embedding)

#-------------------------------------------------------------------------------------------------------------------


def define_param_groups(model, weight_decay, optimizer_name):
    def exclude_from_wd_and_adaptation(name):
        if 'bn' in name:
            return True
        if optimizer_name == 'lars' and 'bias' in name:
            return True

    param_groups = [
        {
            'params': [p for name, p in model.named_parameters() if not exclude_from_wd_and_adaptation(name)],
            'weight_decay': weight_decay,
            'layer_adaptation': True,
        },
        {
            'params': [p for name, p in model.named_parameters() if exclude_from_wd_and_adaptation(name)],
            'weight_decay': 0.,
            'layer_adaptation': False,
        },
    ]
    return param_groups


class SimCLR_pl(pl.LightningModule):
    def __init__(self, config):
        super().__init__()
        self.config = config

        #self.model = AddProjection(config, model=model, mlp_dim=feat_dim)
        self.model = AddProjectionMLP(config)

        self.loss = ContrastiveLoss(config.batch_size, temperature=self.config.temperature)

    def forward(self, X):
        return self.model(X)

    def training_step(self, batch, batch_idx):
        (x1, x2), labels = batch
        z1 = self.model(x1)
        z2 = self.model(x2)
        loss = self.loss(z1, z2)
        self.log('Contrastive loss', loss, on_step=True, on_epoch=True, prog_bar=True, logger=True)
        return loss

    def configure_optimizers(self):
        max_epochs = int(self.config.epochs)
        param_groups = define_param_groups(self.model, self.config.weight_decay, 'adam')
        lr = self.config.lr
        optimizer = Adam(param_groups, lr=lr, weight_decay=self.config.weight_decay)

        print(f'Optimizer Adam, '
              f'Learning Rate {lr}, '
              f'Effective batch size {self.config.batch_size * self.config.gradient_accumulation_steps}')

        scheduler_warmup = LinearWarmupCosineAnnealingLR(optimizer, warmup_epochs=10, max_epochs=max_epochs,
                                                         warmup_start_lr=0.0)

        return [optimizer], [scheduler_warmup]


#-------------------------------------------------------------------------------------------------------------------
"""## Hyperparameters, and configuration stuff"""

# a lazy way to pass the config file
class Hparams:
    def __init__(self):
        self.input_dim = 5000 #number of genes
        self.hidden_dims = [2048, 1024, 512]
        self.embedding_size = 128 
        
        self.epochs = 5 # number of training epochs 300 before changing to 15
        self.seed = 77777 # randomness seed
        self.cuda = True # use nvidia gpu
        # self.img_size = 96 #image shape
        self.save = "./saved_models/" # save checkpoint
        self.load = False # load pretrained checkpoint
        self.gradient_accumulation_steps = 5 # gradient accumulation steps
        self.batch_size = 10 # should be 2*10 in each batch. 
        self.lr = 3e-4 # for ADAm only
        self.weight_decay = 1e-6
        self.temperature = 0.5 # 0.1 or 0.5
        self.checkpoint_path = './scContrastiveLearn.ckpt' # replace checkpoint path here
#-------------------------------------------------------------------------------------------------------------------


"""## Pretraining main logic"""

from pytorch_lightning import Trainer
import os
from pytorch_lightning.callbacks import GradientAccumulationScheduler
from pytorch_lightning.callbacks import ModelCheckpoint
# from torchvision.models import  resnet18

import DataLoader_tensor_sparse as dl
import SCDataset as ds


available_gpus = len([torch.cuda.device(i) for i in range(torch.cuda.device_count())])
save_model_path = os.path.join(os.getcwd(), "saved_models/")
print('available_gpus:',available_gpus)
filename='scContrastiveLearn_adam_'
resume_from_checkpoint = False
train_config = Hparams()

reproducibility(train_config)
save_name = filename + '.ckpt'
#-------------------------------------------------------------------------------------------------------------------
# model = SimCLR_pl(train_config, model=resnet18(pretrained=False), feat_dim=512)
model = SimCLR_pl(train_config)
#-------------------------------------------------DataLoading-------------------------------------------------------
# transform = Augment(train_config.img_size)
# data_loader = get_stl_dataloader(train_config.batch_size, transform)


import anndata as ad
import numpy as np
import scipy
import pandas as pd

normed_counts = "/home/users/syang71/Dataset/Larry_Dataset_normalized/stateFate_inVitro_normed_counts.mtx.gz"  
gene_names = "/home/users/syang71/Dataset/Larry_Dataset_normalized/stateFate_inVitro_gene_names.txt.gz" 
clone_matrix = "/home/users/syang71/Dataset/Larry_Dataset_normalized/stateFate_inVitro_clone_matrix.mtx.gz" 
metadata = "/home/users/syang71/Dataset/Larry_Dataset_normalized/stateFate_inVitro_metadata.txt.gz" 

# load data
normed_counts_mat = scipy.io.mmread(normed_counts).tocsr()
genes = pd.read_csv(gene_names, sep='\t',header=None).to_numpy().flatten()
clone_mat = scipy.io.mmread(clone_matrix).tocsr()
meta_df = pd.read_csv(metadata, sep='\t')

# create full adata
adata = ad.AnnData(normed_counts_mat, obs=meta_df, var=pd.DataFrame(index=genes), dtype=np.float32)
# optimize dtypes
adata.obs['Library'] = adata.obs['Library'].astype('category')
adata.obs['Time point'] = adata.obs['Time point'].astype(int)
adata.obs['Starting population'] = adata.obs['Starting population'].astype('category')
adata.obs['Cell type annotation'] = adata.obs['Cell type annotation'].astype('category')
adata.obs['Well'] = adata.obs['Well'].astype(int)
# assign clone_id
adata.obs['clone_id'] = (clone_mat @ np.arange(1,1+clone_mat.shape[1])) - 1
print("number of lineages: ", len(adata.obs['clone_id'].unique()))
# input data
count_matrix = adata.X
cell_lineage = adata.obs['clone_id'].values.reshape(-1, 1)
count_matrix.shape, cell_lineage.shape
# step 1 generate designed batches
# batchsize = 10
DLoader = dl.SClineage_DataLoader(count_matrix,cell_lineage,batch_size= train_config.batch_size, seed=7)
batch_all, num_batch = DLoader.batch_generator()
# step 2 generate real dataloader
sc_dataset = ds.SCDataset(batches=batch_all)

print("number of batches: ", num_batch)

data_loader = torch.utils.data.DataLoader(dataset=sc_dataset, batch_size=train_config.batch_size, shuffle=False, num_workers=cpu_count()//2)

#-------------------------------------------------------------------------------------------------------------------

accumulator = GradientAccumulationScheduler(scheduling={0: train_config.gradient_accumulation_steps})
checkpoint_callback = ModelCheckpoint(filename=filename, dirpath=save_model_path,
                                        save_last=True, save_top_k=2,monitor='Contrastive loss_epoch',mode='min')

if resume_from_checkpoint:
  trainer = Trainer(callbacks=[accumulator, checkpoint_callback],
                  gpus=available_gpus,
                  max_epochs=train_config.epochs,
                  resume_from_checkpoint=train_config.checkpoint_path)
else:
  trainer = Trainer(callbacks=[accumulator, checkpoint_callback],
                  gpus=available_gpus,
                  max_epochs=train_config.epochs)


trainer.fit(model, data_loader)


trainer.save_checkpoint(save_name)

#-------------------------------------------------------------------------------------------------------------------
"""## Save only backbone weights from Resnet18 that are only necessary for fine tuning"""

# model_pl = SimCLR_pl(train_config, model=resnet18(pretrained=False))
# model_pl = weights_update(model_pl, "SimCLR_ResNet18_adam_.ckpt")

# resnet18_backbone_weights = model_pl.model.backbone
# print(resnet18_backbone_weights)
# torch.save({
#             'model_state_dict': resnet18_backbone_weights.state_dict(),
#             }, 'resnet18_backbone_weights.ckpt')

model_pl = SimCLR_pl(train_config)
model_pl = weights_update(model_pl, "scContrastiveLearn_adam_.ckpt")

baseencoder_backbone_weights = model_pl.model.backbone
print(baseencoder_backbone_weights)
torch.save({
            'model_state_dict': baseencoder_backbone_weights.state_dict(),
            }, 'baseencoder_backbone_weights.ckpt')

#-------------------------------------------------------------------------------------------------------------------

