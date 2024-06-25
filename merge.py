import anndata as ad

def merge(l, output_dir):
    """
    Merge multiple AnnData objects.

    Input:
        l (list of str): List of AnnData object file names to merge.
        output_dir (str): The output directory
    
    Output:
        Merged AnnData object saved as object_all-the-datasets.h5ad in the Objects/ directory.
    """
    # Construct the name for the merged object by concatenating parts of the individual object names
    l_mid = l[1:-1]  # All but the first and last elements in the list
    nw_l = ''
    for i in range(len(l_mid)): 
        nw_l += l_mid[i].replace('.h5ad', '').replace('object', '')

    # Create the final object name by combining the first, middle, and last parts
    obj_name = l[0].replace('.h5ad', '') + nw_l + l[-1].replace('object', '')

    # Initialize dictionaries and lists to store AnnData objects
    d_adata = {}
    l_adata = []
    
    # Read each AnnData object file and store it in the dictionary and list
    for i in range(len(l)):
        # Read the AnnData object from file
        d_adata[f'adata{i}'] = ad.read_h5ad(output_dir + '/Objects/' + l[i])
        l_adata.append(d_adata[f'adata{i}'])

    # Merge all AnnData objects using an outer join to include all genes
    adata_merge = ad.concat(l_adata, join = 'outer')

    # Uncomment the following lines if you want to create a violin plot of the merged object
    # sc.pl.violin(adata_merge, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter = 0.4, size = 1.6, multi_panel = True)

    # Save the merged AnnData object to the specified directory
    adata_merge.write(output_dir + '/Objects/' + obj_name)
