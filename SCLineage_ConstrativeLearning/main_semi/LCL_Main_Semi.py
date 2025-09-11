# General packages
import time
import random
import os
import numpy as np
import anndata as ad
import multiprocessing

# Deep learning packages
import torch
from torch.optim import AdamW
import pytorch_lightning as pl
from pl_bolts.optimizers.lr_scheduler import LinearWarmupCosineAnnealingLR
from torch.utils.data import DataLoader, TensorDataset
from pytorch_lightning.callbacks import Callback, ModelCheckpoint, GradientAccumulationScheduler
from pytorch_lightning.callbacks import EarlyStopping
from pytorch_lightning.loggers import TensorBoardLogger

# Self-written packages
from LCL_Model_Semi import BaseEncoder_ProjHead_MLP, ContrastiveLoss
import General_Dataloader as GD
import SCDataset as ds

# Argument parsing
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Run the contrastive learning model on provided single-cell data.")
    parser.add_argument('--inputFilePath', type=str, help='Anndata for running the algorithm')
    parser.add_argument('--unlabeled_per_batch', type=int, default=5, help='number of unlabeled cells per batch')
    parser.add_argument('--lambda_penalty', type=float, default=1, help='lambda_penalty for semi-supervised learning')
    parser.add_argument('--batch_size', type=int, default=25, help='Batch size for training and validation')
    parser.add_argument('--size_factor', type=float, default=0.3, help='Size factor range from 0 to 1')
    parser.add_argument('--temperature', type=float, default=0.5, help='Temperature parameter for contrastive loss')
    parser.add_argument('--patience', type=int, default=10, help='Number of epochs with no improvement after which training will be stopped')
    parser.add_argument('--min_delta', type=float, default=0.001, help='Minimum change to qualify as an improvement')
    parser.add_argument('--max_epoch', type=int, default=220, help='Maximum number of epochs')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save outputs')
    parser.add_argument('--train_test', default=0, type=int, help='1: split the data for train and validation; 0: otherwise')
    parser.add_argument('--hidden_dims', default=[1024, 256, 64], type=lambda s: [int(item) for item in s.split(',')], help='dimensions of each layer of base encoder. example input: 1024,256,64')
    parser.add_argument('--embedding_size', default=32, type=int, help='the output dimension of projection head')
    parser.add_argument('--resume_from_checkpoint', type=str, help='Path to a checkpoint to resume from', default=None)
    return parser.parse_args()

# Helper functions
def default(val, def_val):
    return def_val if val is None else val

def device_as(t1, t2):
    """Moves t1 to the device of t2"""
    return t1.to(t2.device)

def reproducibility(seed, use_cuda):
    torch.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    np.random.seed(seed)
    if use_cuda:
        torch.cuda.manual_seed(seed)

def load_checkpoint(model, checkpoint_path):
    checkpoint = torch.load(checkpoint_path)
    model.load_state_dict(checkpoint['state_dict'], strict=False)
    print(f'Checkpoint {checkpoint_path} was loaded')
    return model

def define_param_groups(model, weight_decay):
    def exclude_from_wd_and_adaptation(name):
        return 'bn' in name or 'bias' in name
    param_groups = [
        {'params': [p for name, p in model.named_parameters() if not exclude_from_wd_and_adaptation(name)], 'weight_decay': weight_decay},
        {'params': [p for name, p in model.named_parameters() if exclude_from_wd_and_adaptation(name)], 'weight_decay': 0.}
    ]
    return param_groups

# Model definition
class scContraLearn(pl.LightningModule):
    def __init__(self, config):
        super().__init__()
        self.config = config
        # base encoder + projection head
        self.model = BaseEncoder_ProjHead_MLP(config)
        # semi‑supervised loss with NT-Xent + entropy penalty
        self.loss_fn = ContrastiveLoss(
            batch_size=config.batch_size,
            temperature=config.temperature,
            lambda_penalty=config.lambda_penalty
        )

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        # batch = (x1, x2, x_unl) each with an extra leading dim of 1
        x1, x2, x_unl = batch
        x1 = x1.squeeze(0)      # → (N_pos, G)
        x2 = x2.squeeze(0)      # → (N_pos, G)
        x_unl = x_unl.squeeze(0)  # → (N_unl, G) or (0, G)

        # encode + project
        z1 = self.model(x1)   # (N_pos, D)
        z2 = self.model(x2)   # (N_pos, D)
        if x_unl.size(0) > 0:
            z_unl = self.model(x_unl)  # (N_unl, D)
        else:
            z_unl = None

        # compute combined loss
        loss = self.loss_fn(z1, z2, z_unl)

        # logging
        self.log('train_loss',     loss,            on_step=True, on_epoch=True, prog_bar=True)
        self.log('contrast_loss',  self.loss_fn.last_contrastive, on_epoch=True)
        self.log('entropy_penalty', self.loss_fn.last_penalty,     on_epoch=True)
        return loss

    # def validation_step(self, batch, batch_idx):
    #     # only use the labeled pairs in validation
    #     x1, x2, _ = batch
    #     x1, x2 = x1.squeeze(0), x2.squeeze(0)
    #     z1, z2 = self.model(x1), self.model(x2)

    #     # pure contrastive loss (no penalty)
    #     val_loss = ContrastiveLoss(
    #         batch_size=self.config.batch_size,
    #         temperature=self.config.temperature,
    #         lambda_penalty=0.0
    #     )(z1, z2, None)

    #     self.log('val_loss', val_loss, on_epoch=True, prog_bar=True)
    #     return val_loss

    def validation_step(self, batch, batch_idx):
        x1, x2, _ = batch
        x1, x2 = x1.squeeze(0), x2.squeeze(0)
        z1, z2 = self.model(x1), self.model(x2)

        val_loss_fn = ContrastiveLoss(
            batch_size=self.config.batch_size,
            temperature=self.config.temperature,
            lambda_penalty=0.0
        ).to(z1.device)                # <--- move mask buffer to GPU

        val_loss = val_loss_fn(z1, z2, None)
        self.log('val_loss', val_loss, on_epoch=True, prog_bar=True)
        return val_loss

    def validation_epoch_end(self, outputs):
        avg_loss = torch.stack(outputs).mean()
        self.log('avg_val_loss', avg_loss, on_epoch=True, prog_bar=True)

    def configure_optimizers(self):
        # set up AdamW + warmup+cosine schedule
        param_groups = define_param_groups(self.model, self.config.weight_decay)
        optimizer = AdamW(param_groups, lr=self.config.lr, weight_decay=self.config.weight_decay)
        scheduler = LinearWarmupCosineAnnealingLR(
            optimizer,
            warmup_epochs=10,
            max_epochs=self.config.epochs
        )
        return [optimizer], [scheduler]


class LossCallback(Callback):
    def __init__(self):
        super().__init__()
        self.train_losses = []
        self.val_losses = []
        self.avg_val_losses = []

    def on_train_epoch_end(self, trainer, pl_module, unused=None):
        train_loss = trainer.callback_metrics.get('train_loss')
        if train_loss is not None:
            self.train_losses.append(train_loss.item())

    def on_validation_epoch_end(self, trainer, pl_module):
        val_loss = trainer.callback_metrics.get('val_loss')
        avg_val_loss = trainer.callback_metrics.get('avg_val_loss')
        if val_loss is not None:
            self.val_losses.append(val_loss.item())
        if avg_val_loss is not None:
            self.avg_val_losses.append(avg_val_loss.item())

class SaveCheckpointCallback(Callback):
    def on_train_epoch_end(self, trainer, pl_module):
        save_path = os.path.join(trainer.checkpoint_callback.dirpath, "scContrastiveLearn_last.ckpt")
        trainer.save_checkpoint(save_path)
        # print(f"Checkpoint saved at {save_path}")


# Configuration class
class Hparams:
    def __init__(self, args):
        self.input_dim = 2000  # number of genes
        # self.hidden_dims = [1024, 256, 64]
        # self.embedding_size = 32  # size of the output embeddings
        self.hidden_dims = args.hidden_dims
        self.embedding_size = args.embedding_size # size of the output embeddings
        
        self.unlabeled_per_batch = args.unlabeled_per_batch
        self.lambda_penalty = args.lambda_penalty
        self.epochs = args.max_epoch  # number of training epochs
        self.batch_size = args.batch_size
        self.size_factor = args.size_factor
        self.temperature = args.temperature
        self.train_test_ratio = 0.8  # ratio for splitting train/test data
        self.patience = args.patience
        self.min_delta = args.min_delta

        self.seed = 3407  # seed for numpy, PyTorch (reproducibility)
        self.batch_seed = 17  # seed for creating batches
        self.random_seed = 42  # seed for the Python 'random' module (data shuffling)
        self.train_test_seed = 42  # seed for splitting train/test data

        self.gradient_accumulation_steps = 5  # number of gradient accumulation steps
        self.lr = 3e-4  # learning rate
        self.weight_decay = 1e-6  # weight decay for regularization
        self.cuda = True  # whether to use CUDA
        self.load = False  # whether to load a saved model
        self.train_test = bool(args.train_test)
        self.save = args.output_dir+"/saved_models/"  # path to save models
        self.checkpoint_path = args.output_dir+'/scContrastiveLearn.ckpt'  # path to checkpoint file
        self.out_dir = args.output_dir
        self.file_path = args.inputFilePath


def prepare_data_loaders(config):
    # 1) pull in labeled + unlabeled batches
    batch_all_label, batch_all_unlabel, num_batches, lineage_array_label, lineage_array_unlabel = \
        GD.General_DataLoader(
            config.file_path,
            config.batch_size,
            config.size_factor,
            config.unlabeled_per_batch,
            config.batch_seed
        )

    # Determine the number of available CPU cores
    num_workers = min(12, multiprocessing.cpu_count())
    print("num_workers(number of available CPU cores):", num_workers)
    print(f"Total batches: {num_batches}")
    print(f"Pos‑pairs per batch: {config.batch_size}, unlabeled per batch: {config.unlabeled_per_batch}")

    # 2) train/val split at the batch‑level

    np.save(config.out_dir + f'/all_lineage_label_bs{config.batch_size}_tau{config.temperature}.npy',lineage_array_label)
    all_unlab = np.concatenate(list(lineage_array_unlabel.values()), axis=0)
    np.save(config.out_dir + f'/all_lineage_unlabel_bs{config.batch_size}_tau{config.temperature}.npy',all_unlab)


    all_keys = list(batch_all_label.keys())
    if config.train_test:
        n_train = int(config.train_test_ratio * num_batches)
        train_keys = random.sample(all_keys, n_train)
        val_keys   = [k for k in all_keys if k not in train_keys]

        # build per‑split dicts
        train_label   = {k: batch_all_label[k]   for k in train_keys}
        train_unlabel = {k: batch_all_unlabel[k] for k in train_keys}
        val_label     = {k: batch_all_label[k]   for k in val_keys}
        val_unlabel   = {k: batch_all_unlabel[k] for k in val_keys}

        print("Training/Validation split:")
        print(f" number of train batches: {len(train_keys)}, number of val batches: {len(val_keys)}")

        train_ds = ds.SCDataset(train_label, train_unlabel)
        val_ds   = ds.SCDataset(val_label,   val_unlabel)

        return {
            'train': DataLoader(train_ds, batch_size=1, shuffle=False,  num_workers=num_workers),
            'val':   DataLoader(val_ds,   batch_size=1, shuffle=False, num_workers=num_workers)
        }

    else:
        # use all data
        all_ds = ds.SCDataset(batch_all_label, batch_all_unlabel)
        return {
            'all': DataLoader(all_ds, batch_size=1, shuffle=False, num_workers=num_workers)
        }

def save_model_and_losses(trainer, loss_callback, config):
    save_name = os.path.join(config.save, 'scContrastiveLearn_.ckpt')
    trainer.save_checkpoint(save_name)

    # Save training losses
    np.save(os.path.join(config.out_dir, f'train_losses_bs{config.batch_size}_tau{config.temperature}.npy'), np.array(loss_callback.train_losses))

    if config.train_test:
        # Save validation losses if train_test is True
        np.save(os.path.join(config.out_dir, f'val_losses_bs{config.batch_size}_tau{config.temperature}.npy'), np.array(loss_callback.val_losses))
        np.save(os.path.join(config.out_dir, f'avg_val_losses_bs{config.batch_size}_tau{config.temperature}.npy'), np.array(loss_callback.avg_val_losses))



def main():
    args = get_args()
    train_config = Hparams(args)
    # Ensure directory exists
    if not os.path.exists(train_config.save):
        os.makedirs(train_config.save)

    reproducibility(train_config.seed, train_config.cuda)

    print("-------------------------------INFO-------------------------------")
    print("Anndata Info: ", train_config.file_path)
    print("batch_size: ", train_config.batch_size)
    print("size_factor: ", train_config.size_factor)
    print("temperature: " , train_config.temperature)
    print("number of epochs: ", train_config.epochs)
    print("train_test_ratio: ", train_config.train_test_ratio)
    print("input_dim: ", train_config.input_dim)
    print("hidden_dims: ", train_config.hidden_dims)
    print("embedding_size: ", train_config.embedding_size)

    model = scContraLearn(train_config)
    data_loaders = prepare_data_loaders(train_config)

    # call back
    loss_callback = LossCallback()
    
    checkpoint_callback = ModelCheckpoint(
        dirpath=train_config.save,
        filename='scContrastiveLearn_{epoch:02d}',
        save_last=True,
        save_top_k=2,
        monitor='avg_val_loss' if train_config.train_test else 'train_loss',
        mode='min'
    )

    early_stopping_callback = EarlyStopping(
        monitor='avg_val_loss' if train_config.train_test else 'train_loss',
        patience=train_config.patience,
        min_delta=train_config.min_delta,
        mode='min',
        verbose=True
    )

    tb_logger = TensorBoardLogger(save_dir=train_config.out_dir,name="lightning_logs",version=None)

    trainer = pl.Trainer(
        logger=tb_logger,
        callbacks=[GradientAccumulationScheduler({0: train_config.gradient_accumulation_steps}), checkpoint_callback, loss_callback, early_stopping_callback, SaveCheckpointCallback()],
        gpus=torch.cuda.device_count(),
        max_epochs=train_config.epochs,
        check_val_every_n_epoch=1 if train_config.train_test else None,
        resume_from_checkpoint=args.resume_from_checkpoint
    )


    if train_config.train_test:
        trainer.fit(model, train_dataloaders=data_loaders['train'], val_dataloaders=data_loaders['val'])
    else:
        trainer.fit(model, train_dataloaders=data_loaders['all'])

    save_model_and_losses(trainer, loss_callback, train_config)

    print("trainning done")
    print("")
    print("-------------------------------Feature Extractin-------------------------------")
    ##-----------------------Extract Features generated by the base encoder for cell pairs-----------------------------
    model.eval()  # Set the model to evaluation mode
    model.to('cuda' if torch.cuda.is_available() else 'cpu')  # Move model to the appropriate device

    ##------------------------------Extract Features generated by the base encoder for cells----------------------------
    print("#--------------------------------------------Feature Extracting(pairs)------------------------------------------------------")
    model.eval()  # Set the model to evaluation mode
    model.to('cuda' if torch.cuda.is_available() else 'cpu')  # Move model to the appropriate device


    adata_subset = ad.read_h5ad(train_config.file_path)
    count_matrix = adata_subset.X
    try:
        count_matrix_arr = count_matrix.toarray()
    except AttributeError:    
        count_matrix_arr = count_matrix

    count_matrix_th = torch.from_numpy(count_matrix_arr)
    dataset_cell = TensorDataset(count_matrix_th)
    data_loader_all = torch.utils.data.DataLoader(dataset_cell, batch_size=train_config.batch_size, shuffle=False, num_workers=1,drop_last=False)
    print("num of batches for all cells (not cell pairs):", len(data_loader_all))


    features_list_Z = []

    for batch in data_loader_all:  
        Z = batch[0]
        Z = Z.to(next(model.parameters()).device)  
        with torch.no_grad():  # No need to compute gradients
            batch_features_Z = model.model.get_features(Z)  # Extract features
            
        features_list_Z.append(batch_features_Z.cpu().detach().numpy())  # Store features as NumPy array

    # Concatenate all batch features into a single NumPy array
    features_Z = np.concatenate(features_list_Z, axis=0)


    print("Shape of the feature representation generated by the base encoder:", features_Z.shape)
    np.save(train_config.out_dir+ f'/scBaseEncoderFeat_Z_bs{train_config.batch_size}_tau{train_config.temperature}.npy', features_Z)


if __name__ == "__main__":
    start_time = time.time()
    print("start time:", start_time)
    main()
    end_time = time.time()
    print("end time:", end_time)
    print(f"Execution time: {round((end_time - start_time)/3600,2)} hours")