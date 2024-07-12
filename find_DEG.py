import os
import scanpy as sc
import anndata as ad
import re

def find_DEG(obj, dir):
    """
    Create UMAPs from an object and save them. Need a list of genes which user have to create in dir Analysis with one gene per line.

    Input:
        obj (str): Object selected by the user
        dir (str): Processing dir

    Output: UMAPs in Analysis/name_dir_of_objects/
            It creates 2 UMAPs : one for datasets and louvain, the second with list of genes
            Generate a log file (umap_parameters.tab) containing information on the UMAP (parameters ...)
    """
    # Create folder
    if not os.path.isdir('Analysis/' + dir + '/DEG'):
        os.mkdir('Analysis/' + dir + '/DEG')

    # Read the AnnData object
    adata = ad.read_h5ad('Processing/' + dir + '/Objects/Objects_clusters/' + obj)
    name = obj.replace('.h5ad', '').replace('o', 'O')

    # Create folder
    if not os.path.isdir('Analysis/' + dir + f'/DEG/{name}'):
        os.mkdir('Analysis/' + dir + f'/DEG/{name}')

    # Take the number of clusters
    l = adata.obs['louvain'].cat.categories.tolist()
    sc.tl.rank_genes_groups(adata, 'louvain')

    # Save the DEG
    for cluster in l:
        df = sc.get.rank_genes_groups_df(adata, group=str(cluster))
        filtered_df = df[df['pvals_adj'] < 0.05]
        filtered_df.to_csv('Analysis/' + dir + f'/DEG/{name}/cluster_{cluster}.tab', sep = '\t')

    print('Done')