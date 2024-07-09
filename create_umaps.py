import os
import scanpy as sc
import anndata as ad
import matplotlib as mpl
from matplotlib import pyplot as plt
from datetime import datetime
import tkinter as tk
import tkfilebrowser
from tkinter import simpledialog, messagebox

# Define analysis parameters here if you want to change
resolution = 0.8  # resolution for clustering, more high = more clusters, low = less clusters
n_neighbors = 50  # number of neighbors for compute the neighborhood graph
n_pcs = 20  # number of principal components analysis
show_pca = False  # open a window and show the pca plot WARNING while this window is open the script is in pause
show_elbow_plot = False  # same as show pca but show elbow plot
max_iter_harmony = 10  # maximum iteration of harmony

def create_umaps(obj, dir):
    """
    Create UMAPs from an object and save them. Need a list of genes which user have to create in dir Analysis with one gene per line.

    Input:
        obj (str): Object selected by the user
        dir (str): Processing dir

    Output: UMAPs in Analysis/name_dir_of_objects/
            It creates 2 UMAPs : one for datasets and louvain, the second with list of genes
            Generate a log file (umap_parameters.tab) containing information on the UMAP (parameters ...)
    """
    sc.settings.set_figure_params(dpi=150, facecolor='white', fontsize=12)
    
    # Ensure the 'Analysis' directory and the specific analysis directory exist
    if not os.path.isdir('Analysis'):
        os.mkdir('Analysis')
    if not os.path.isdir('Analysis/' + dir):
        os.mkdir('Analysis/' + dir)

    # Tkinter root window for input dialog
    root = tk.Tk()
    root.withdraw()

    # User inputs for the UMAP parameters
    theta = simpledialog.askstring('Input', 'Enter the theta (None/1/2/...):')
    # Convert theta to appropriate type to avoid errors with harmony
    if theta == 'None':
        theta = None
    else:
        theta = int(theta)

    nclust = simpledialog.askstring('Input', 'Enter the number of clusters (None/1/2/3...):')

    # Convert nclust to appropriate type to avoid errors with harmony
    if nclust == 'None':
        nclust = None
    else:
        nclust = int(nclust)

    messagebox.showinfo('Information', 'Select the file which contains list gene (1 gene per line)')
    file_l_gene = tkfilebrowser.askopenfilename(initialdir = 'Analysis/', title = 'Select the list gene file')

    # Read the list of genes from the specified file
    l_gene = []
    with open(file_l_gene, 'r') as f:
        for lig in f:
            lig = lig.rstrip()
            l_gene.append(lig)
    
    # Read the AnnData object
    adata = ad.read_h5ad('Processing/' + dir + '/Objects/' + obj)

    # Normalize the data
    sc.pp.normalize_total(adata, target_sum = 1e4)

    # Log-transform the data
    sc.pp.log1p(adata)

    # Store the raw data
    adata.raw = adata

    # Scale the data
    sc.pp.scale(adata)

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata)

    # Run PCA
    sc.tl.pca(adata)

    # Plot PCA if requested
    if show_pca:
        sc.pl.pca(adata, size = 18, na_color = 'black')

    # Plot elbow plot if requested
    if show_elbow_plot:
        sc.pl.pca_variance_ratio(adata)

    # Integrate data using Harmony
    sc.external.pp.harmony_integrate(adata, 'dataset', theta = theta, max_iter_harmony = max_iter_harmony, nclust = nclust)

    # Save the Harmony-corrected PCA results
    adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

    # Compute the neighborhood graph
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, n_pcs = n_pcs)

    # Cluster the cells using the Louvain algorithm
    sc.tl.louvain(adata, resolution = resolution)

    # Compute UMAP
    sc.tl.umap(adata)

    # Ensure traceability by creating a log file if it doesn't exist
    if not os.path.exists('Analysis/' + dir + '/umaps_parameters.tab'):
        with open('Analysis/' + dir + '/umaps_parameters.tab', 'a') as f:
            # Get the current date and time
            date = datetime.now()
            date_str = date.strftime('%Y-%m-%d %H:%M:%S')
            # Write the header to the log file
            f.write(date_str + '\n' + 'umap_number' + '\t' + 'object' + '\t' 'theta' + '\t' + 'nb_clust' + '\t' + 'resolution' + 
                    '\t' + 'nb_neighbors' + '\t' + 'n_pcs' + '\n')

    # Generate a unique filename for the UMAP datasets Louvain plot
    i = 1
    name_umap_datasets_louvain = 'umap_datasets_louvain_1.pdf'
    while os.path.exists('Analysis/' + dir + '/' + name_umap_datasets_louvain):
        i += 1
        name_umap_datasets_louvain = f'umap_datasets_louvain_{i}.pdf'

    # Plot and save the UMAP datasets Louvain plot
    with plt.rc_context():
        sc.pl.umap(adata, color=['dataset', 'louvain'], show=False)
        plt.savefig('Analysis/' + dir + '/' + name_umap_datasets_louvain)

    # Generate a unique filename for the UMAP genes plot
    i = 1
    name_umap_genes = 'umap_genes_1.pdf'
    while os.path.exists('Analysis/' + dir + '/' + name_umap_genes):
        i += 1
        name_umap_genes = f'umap_genes_{i}.pdf'

    # Plot and save the UMAP genes plot
    with plt.rc_context():
        sc.pl.umap(adata, color=l_gene, color_map=mpl.cm.Reds, show=False)
        plt.savefig('Analysis/' + dir + '/' + name_umap_genes)

    # Append filter information to the log file
    with open('Analysis/' + dir + '/umaps_parameters.tab', 'a') as f:
        f.write(f'umap_{i}' + '\t' + dir + '/Objects/' + obj + '\t' + str(theta) + 
                '\t' + str(nclust) + '\t' + str(resolution) + '\t' + str(n_neighbors) + '\t' + str(n_pcs) + '\n')
