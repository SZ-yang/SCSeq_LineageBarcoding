import anndata as ad
import numpy as np
import scipy
import pandas as pd
import scanpy as sc
import SCSeq_LineageBarcoding2.SCSeq_LineageBarcoding.SCLineage_ConstrativeLearning.main.DataLoader_tensor_sparse as dl
import SCDataset as ds

def Larry_DataLoader(num_genes, batch_size, batch_seed):

#---------------------------------------------------------Load the matrixs-----------------------------------------------------

    normed_counts = "/Users/apple/Desktop/KB/Dataset1/stateFate_inVitro_normed_counts.mtx.gz"  #snakemake.input['normed_counts']
    gene_names = "/Users/apple/Desktop/KB/Dataset1/stateFate_inVitro_gene_names.txt.gz" #snakemake.input['gene_names']
    clone_matrix = "/Users/apple/Desktop/KB/Dataset1/stateFate_inVitro_clone_matrix.mtx.gz" #snakemake.input['clone_matrix']
    metadata = "/Users/apple/Desktop/KB/Dataset1/stateFate_inVitro_metadata.txt.gz" #snakemake.input['metadata'] 

    # load data
    normed_counts_mat = scipy.io.mmread(normed_counts).tocsr()
    genes = pd.read_csv(gene_names, sep='\t',header=None).to_numpy().flatten()
    clone_mat = scipy.io.mmread(clone_matrix).tocsr()
    meta_df = pd.read_csv(metadata, sep='\t')



#-------------------------------Get num_genes of highly expressed genes from the orginal data(all cells)--------------------
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

    # get 2000 genes from the 130887(all) cells
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,n_top_genes=num_genes)
    

#---------------------------------Creat the subset adata with trimmed number of genes and number of cells--------------------
    # create full adata
    adata_cp = ad.AnnData(normed_counts_mat, obs=meta_df, var=pd.DataFrame(index=genes), dtype=np.float32)

    # optimize dtypes
    adata_cp.obs['Library'] = adata_cp.obs['Library'].astype('category')
    adata_cp.obs['Time point'] = adata_cp.obs['Time point'].astype(int)
    adata_cp.obs['Starting population'] = adata_cp.obs['Starting population'].astype('category')
    adata_cp.obs['Cell type annotation'] = adata_cp.obs['Cell type annotation'].astype('category')
    adata_cp.obs['Well'] = adata_cp.obs['Well'].astype(int)
    # assign clone_id
    adata_cp.obs['clone_id'] = (clone_mat @ np.arange(1,1+clone_mat.shape[1])) - 1

    # remove the cells that belong to the clone_id with few cells
    value_counts = adata.obs['clone_id'].value_counts()
    frequency_dict = {}
    for value, count in value_counts.items():
        if count in frequency_dict:
            frequency_dict[count].append(value)
        else:
            frequency_dict[count] = [value]

    clone_for_remove = frequency_dict[81585]+frequency_dict[2]+frequency_dict[3]+frequency_dict[4] 
    adata_subset = adata_cp[~adata_cp.obs['clone_id'].isin(clone_for_remove)]
    print("adata_subset.obs.shape:", adata_subset.obs.shape)

    # trim the number of genes to the required number using the highly variable gene calculated from the original adata (all cells)
    hvgene = (adata.var.highly_variable[adata.var.highly_variable==True]).index
    print("number of the highly variable genes:", len(hvgene))
    adata_subset = adata_subset[:,hvgene]
    print("adata_subset.X.shape:", adata_subset.X.shape)

#----------------------------------------------------Generate Batches--------------------------------------------------------
    # input data
    count_matrix = adata_subset.X
    cell_lineage = adata_subset.obs['clone_id'].values.reshape(-1, 1)
    # count_matrix.shape, cell_lineage.shape

    # step 1 generate designed batches
    # batchsize = 10
    DLoader = dl.SClineage_DataLoader(count_matrix,cell_lineage,batch_size= batch_size, seed=batch_seed)
    batch_all, num_batch, lineage_info = DLoader.batch_generator()
    # step 2 generate real dataloader
    # sc_dataset = ds.SCDataset(batches=batch_all)

    print("number of batches: ", num_batch)
    print("total number of pairs: ", num_batch*batch_size)

    return batch_all, lineage_info, num_batch

if __name__ == "__main__":
    import torch
    larry_dataset, lineage_info = Larry_DataLoader(2000,10)
    data_loader = torch.utils.data.DataLoader(dataset=larry_dataset, batch_size=10, shuffle=False, num_workers=1) #num_workers=cpu_count()//2

    # save the lineage info for UMAP plotting
    np.save('/home/users/syang71/kzlinlab/projects/scContrastiveLearn/git/scContrastiveLearn_Joshua/lineage_info_322.npy', lineage_info)