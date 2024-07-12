# Single-cell RNA-seq data integration

Internship project to integrate differents datasets of single-cell RNA-seq. These scripts can process data in the form of one cell per file (often Smart-seq) or in format of 10X (3 files : barcodes.tsv, genes.tsv, matrix.mtx). It could take all the technologies of single cell but in this 2 formats.

## Installation and Requirements

You have to download these directory and put your data in a directory named **Data**. For format Smart_seq, named a dir **NumberOfDataset_Author_counts** and for 10X, **NumberOfDataset_Author_10X**. 
>[!IMPORTANT]
>If you don't respect named of Data folder, it will not working.

These scripts use multiple packages. For installing, you can use these command lines in your terminal (on vs code terminal for exemple) or you can use conda (see below).
They also use packages pywin32 (**only if you have Windows**).
Use Python version >= 3.11

```
pip install scanpy==1.10.2
pip install harmonypy==0.0.9
pip install igraph==0.10.13
pip install louvain==0.8.2
pip install tkfilebrowser==2.3.2
pip install pywin32==306
```

You can also use conda thanks to the yaml file.
```
conda env create -n scrnaseqint --file scrnaseqint.yaml
```
>[!NOTE]
>It's possible that when you execute `main.py`, an error like `_tkinter.TclError: Item ... already exists`, you can refere here : https://github.com/hwstar/BOMtools/issues/2


## Contents

- `apply_filters.py`: function to apply filters on anndata object on the minimum of genes per cell, the minimum of cells per genes and the maximum of % mitochondrial genes, generate a violin plot and ask the user if he want to save the filtered object
- `create_anndata_object.py`: function to create anndata object to save data and generate a violin plot
- `create_umaps.py`: function to create UMAPs from a merged object
- `find_DEG.py`: function to find differentially expressed genes per cluster from UMAPs
- `gene_names.tab`: ouput of `extraction_gene_name.sh` (in the folder Appendix scripts) with in column 1 gene ids and column 2 gene names, used by `generation_matrix.py`
- `generation_matrix.py`: function to create counts matrix from datasets with 1 cell per file, with columns represent cells and rows corespond to genes, use `gene_names.tab`
- `main.py`: main script which call function and interact with the user, use `apply_filters.py`, `create_anndata_object.py`, `create_umaps.py`, `find_DEG.py`, `generation_matrix.py`, `main.py`, `merge.py`
- `merge.py`: function to merge multiple anndata objects
- `scrnaseqint.yaml`: file to create conda environment

## Data Processing Script

This script interacts with the user to process and analyze single-cell RNA sequencing datasets (Smart-seq2 or 10X). Follow the process to generate matrices, create objects, filter data, merge objects, create UMAPs and find differentially expressed genes per cluster from UMAPs.

## Usage

You have to put your data in a folder named `Data`, and named the folder X_datasetname_10X (for 10X format) or X_datasetname_counts (for matrix format) with X a number. 
If you use VS code, you can Open Folder to be in the right folder. Execute `main.py` script and follow the processus to select the desired operation.

```
python main.py
```

## Main Menu Options

At the beginning, a window appear and show folders in Processing. Here, you can create a folder or select an existing folder to save processing files for a better tracability. (More details in the tutorial folder)

1. **Smart-seq2**
2. **10X**
3. **Apply filters**
4. **Merge objects**
5. **Create UMAPs**
6. **Find differentially expressed genes per cluster from UMAPs**

### Smart-seq2 format objects creation

1. **Generation of matrix in the case of a file per cell**
    - Displays the list of folders containing counts (e.g., `1_Deng_counts` in the Data folder). This applies to datasets with one cell per file in **.tab**.
    - The script executes the `generation_matrix` function from the `generation_matrix.py` script to merge all files of the dataset into a matrix with genes in rows and cells in columns and replacing gene identifiers with gene names. The resulting file will be saved in the `Data` folder (e.g., `matrix_1_Deng.tab`). If matrices already exist, click on Cancel if your files already in matrix format.

2. **Creation of objects (AnnData)**
    - Show folder containing matrix (`Data`)
    - Choose one or more matrix
    - Creates objects in the `Objects` folder (e.g., `object_1_Deng_ori.h5ad`) and generates a violin plot in the `Plots` folder (e.g., `Plots_1_Deng`) using the `create_anndata_object` function from the `create_anndata_object.py` script.

![anndata object](Images/anndata_object.png)

### 10X format objects creation

- Choose one or more folders for process 10X
- Creates objects using the `create_anndata_object` function (similar to Smart-seq2).

### Filtering Data

- Choose one or more objects to add filters
- Requests filter parameters (minimum genes per cell, minimum cells per gene, and maximum percentage of mitochondrial genes).
- Executes the `apply_filters` function from the `apply_filters.py` script.
- Displays the number of cells and genes before and after filtering 
- Asks if the user wants to save the filtered object (e.g. `object_1_Deng_filtered_1.h5ad`) in `Objects/Objects_ori` and saves a violin plot with the filter parameters in the name. Moreover, a log file `filters_appliled.tab` is created with all the parameters choosen to a better tracability.

### Merge Objects

- Displays objects folder
- Choose one or more objects in `Objects/Objects_filtered` to merge (the best is to merge one object by one, not an already merge with another because the script take all genes and some genes can miss if you merged an already merge object with an other moreover it gives errors because of 0)
- Merges objects using the `merge` function from the `merge.py` script and saves the merged object in `Objects/Objects_merged` (e.g. `object_merged_1.h5ad`) and a log file `objects_merged.tab` with the object number and all the object names merged.

### Create UMAPs
- Displays objects folder
- Choose one object from `Objects/Objects_merged` to create UMAPs.
- Choose parameters : theta, number of clusters and a file containing list of genes (one per line)
- Creates 2 UMAPs in `Analysis/folder_name_at_the_beginning/UMAPs` and a log file (`umpas_parameters.tab`) for tracability using `create_umaps.py`.
- Ask if you want to save object file withs clusters to find differentially expressed genes per cluster from UMAPs. If yes, it saves the object `object_clusters_X` (`Objects/Objects_clusters`) with X the number of the UMAP 
- Ask if you want to create other UMAPs
- If you want to change more precise parameters you can directly modify in the script (number of neighbors, number of PCs, show PCA plot, show Elbow plot)

## Find differentially expressed genes per cluster from UMAPs
- Displays objects folder
- Choose one object from `Objects/Objects_clusters`
- Find differentially expressed genes per cluster against others clusters from UMAPs and generates one table per cluster (e.g. `cluster_0.tab`) in `Analysis/folder_name_at_the_beginning/DEG/Object_clusters_X` with X the number of the corresponding UMAPs

## Notes

- A detailled tutorial is available in the folder Tutorial
