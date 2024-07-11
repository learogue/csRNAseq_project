import anndata as ad
import os

def merge(l, output_dir):
    """
    Merge multiple AnnData objects.

    Input:
        l (list of str): List of AnnData object file names to merge.
        output_dir (str): The output directory
    
    Output:
        Merged AnnData object saved as object_all-the-datasets.h5ad in the Objects/ directory.
    """
    # Ensure traceability by creating a log file if it doesn't exist
    if not os.path.exists(output_dir + '/Objects/objects_merged.tab'):
        with open(output_dir + '/Objects/objects_merged.tab', 'a') as f:
            # Write the header to the log file
            f.write('object_number' + '\t' + 'objects_selected' + '\n')

    # Generate a unique filename for the object
    a = 1
    obj_name = 'object_merged_1.h5ad'
    while os.path.exists(output_dir + '/Objects/' + obj_name):
        a += 1
        obj_name = f'object_merged_{a}.h5ad'

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

    # Save the merged AnnData object to the specified directory
    adata_merge.write(output_dir + '/Objects/' + obj_name)

    # Append object information to the log file
    with open(output_dir + '/Objects/objects_merged.tab', 'a') as f:
        f.write(f'object_merged_{a}' + '\t' + ' '.join(l) + '\n')
