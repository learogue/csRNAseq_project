import scanpy as sc
import anndata as ad
from matplotlib import pyplot as plt
from datetime import datetime
import os
import tkinter as tk
from tkinter import simpledialog, messagebox

def apply_filters(file, output_dir):
    """
    Apply filters to AnnData object. User chooses filters on:
        - Minimum number of genes per cell
        - Minimum number of cells per gene
        - Maximum percentage of mitochondrial genes

    Input:
        file (str): AnnData object file name
        output_dir (str): The output directory
    
    Output:
        Violin plot with applied filters saved in Processing/output_dir/Plots/Plots_number-of-dataset_author/
        If user wants to save the filtered data: AnnData object with _filtered_X suffix
        Generate a log file (filters_applied.tab) containing information on the filters applied
    """
    # Create folder
    if not os.path.isdir(output_dir + '/Objects/Objects_filtered'):
        os.mkdir(output_dir + '/Objects/Objects_filtered')    
    
    # Ensure traceability by creating a log file if it doesn't exist
    if not os.path.exists(output_dir + '/Objects/Objects_filtered/filters_applied.tab'):
        with open(output_dir + '/Objects/Objects_filtered/filters_applied.tab', 'a') as f:
            # Get the current date and time
            date = datetime.now()
            date_str = date.strftime('%Y-%m-%d %H:%M:%S')
            # Write the header to the log file
            f.write(date_str + '\n' + 'object' + '\t' + 'nb_cells_before' + '\t' + 'nb_genes_before' + '\t' + 'min_genes_per_cells' + 
                    '\t' + 'min_cells_per_genes' + '\t' + 'max_pct_mt' + '\t' + 'nb_cells_after' + '\t' + 'nb_genes_after' + '\n')
            
    # Extract the name of the dataset from the file name
    name = file.replace('object_', '').replace('_ori.h5ad', '')

    # Read the AnnData object from the file
    adata = ad.read_h5ad(output_dir + '/Objects/Objects_ori/' + file)

    # Print the name of the object and the number of cells and genes before applying filters
    nb_cells_before = adata.n_obs
    nb_genes_before = adata.n_vars

    # User inputs for the filter criteria
    min_genes = simpledialog.askinteger('Input', f'For {name}, Enter a minimum number of genes per cell:')
    min_cells = simpledialog.askinteger('Input', f'For {name}, Enter a minimum number of cells per gene:')
    pct_mt = simpledialog.askinteger('Input', f'For {name}, Enter a maximum percentage of mitochondrial genes:')

    # Apply the filters to the AnnData object
    sc.pp.filter_cells(adata, min_genes = min_genes)
    sc.pp.filter_genes(adata, min_cells = min_cells)
    adata = adata[adata.obs.pct_counts_mt < pct_mt, :].copy()  # Keep only cells with mitochondrial gene percentage below the threshold

    # Number of cells and genes after applying filters
    nb_cells_after = adata.n_obs
    nb_genes_after = adata.n_vars

    # Display results before and after filters in a message box
    result_message = (
        f'Results for {name}:\n'
        f'Number of cells before filters: {nb_cells_before}\n'
        f'Number of genes before filters: {nb_genes_before}\n'
        f'Number of cells after filters: {nb_cells_after}\n'
        f'Number of genes after filters: {nb_genes_after}\n')
    messagebox.showinfo('Filter Results', result_message)

    # Generate and save a violin plot showing the filtered data
    plot_dir = os.path.join(output_dir, 'Plots', f'Plots_{name}')
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    with plt.rc_context(): 
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter = 0.4, size = 1.6, multi_panel = True, show = False)
        plt.savefig(os.path.join(plot_dir, f'violin_plot_min_genes_{min_genes}_min_cells_{min_cells}.pdf'), bbox_inches = 'tight')

    # Ask the user if they want to save the filtered AnnData object
    save = messagebox.askyesno('Save Filtered Object', 'Do you want to save the new object with these filters?')
    if save:
        # Generate a unique filename for the filtered AnnData object
        i = 1
        name_file_output = file.replace('ori.h5ad', 'filtered_1.h5ad')
        while os.path.exists(output_dir + '/Objects/Objects_filtered/' + name_file_output):
            i += 1
            name_file_output = file.replace('ori.h5ad', f'filtered_{i}.h5ad')
        
        # Save the filtered AnnData object
        adata.write(output_dir + '/Objects/Objects_filtered/' + name_file_output)

        # Append filter information to the log file
        with open(output_dir + '/Objects/Objects_filtered/filters_applied.tab', 'a') as f:
            f.write(name_file_output + '\t' + str(nb_cells_before) + '\t' + str(nb_genes_before) + 
                    '\t' + str(min_genes) + '\t' + str(min_cells) + '\t' + str(pct_mt) + '\t' + str(adata.n_obs) + '\t' + str(adata.n_vars) + '\n')

    print('Done')