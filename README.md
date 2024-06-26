# Single-cell RNA-seq datas integration

Internship project to integrate differents datasets of single-cell RNA-seq. These scripts can process datas in the form of one cell per file (often Smart-seq) or in format of 10X (3 files : barcodes.tsv, genes.tsv, matrix.tsv). It could take all the technologies of single cell but in this 2 formats.
I also analyze datasets. I take 6 datasets to perform a standard because of poor number cells. I worked on scRNA-seq from mouse embryos at stage E4.5 (peri-implantation) and I was interested in Lefty1 to analyze cells which differentially expressed this genes which have a important role in the establishement of anterio-posterior axis.

## Installation

You have to download these directory and put your datas in a directory named Data. For format Smart_seq, named a dir **NumberOfDataset_Author_counts** and for 10X, **NumberOfDataset_Author_10X**. See the scheme below which explained the required organisation of the directory.

  ![dir organisation](Images/organisation_scheme.png)

These scripts use the packages scanpy. For installing scanpy, you can use this command line in your terminal (on vs code)

```
  pip install scanpy
```

## Contents

- `apply_filters.py`: function to apply filters on anndata object on the minimum of genes per cell, the minimum of cells per genes and the maximum of % mitochondrial genes, generate a violin plot and ask the user if he want to save the filtered object (_filterd.h5ad)
- `extraction_gene_name.sh`: script used to generate `gene_names.tab` to have gene ids and gene names in order to replace gene ids by names in matrix
- `gene_names.tab`: ouput of `extraction_gene_name.sh` with in column 1 gene ids and column 2 gene names
- `create_anndata_object.py`: function to create anndata object to stock datas, generate a violin plot
- `generation_matrix.py`: function to create counts matrix from datasets with 1 cell per file, with columns represent cells and rows corespond to genes
- `main.py`: main script which call function and interact with the user
- `merge.py`: function to merge multiple anndata object, so it merged datasets, use `apply_filters.py`, `create_anndata_object.py`, `generation_matrix.py`, `main.py`, `merge.py`
- `modify_ids_names.py` : script in case of concatenation of files from Smart-seq2 already performed, and the 1st column correspond to gene IDs, this script change gene ids in gene names and need of `gene_names.tab`

## Data Processing Script

This script interacts with the user to process and merge single-cell RNA sequencing datasets (Smart-seq2 or 10X). Follow the prompts to generate matrices, create objects, filter data, and merge objects.

## Usage

Run the script and follow the prompts to select the desired operation.

### Main Menu Options

1. **Smart-seq2**
2. **10X**
3. **Both Smart-seq2 and 10X**
4. **Merge existing objects only**

### Smart-seq2 (answer 1 or 3)

1. **Generation of matrices in the case of a file per cell**
    - Displays the list of folders containing counts (e.g., `1_Deng_counts` in the Data folder). This applies to datasets with one cell per file in tabular format.
    - The script asks for the number of directories to process.
    - Enter the folder names one by one.
    - The script executes the `generation_matrix` function from the `generation_matrix.py` script to merge all files of the dataset into a matrix with genes in rows and cells in columns, replacing gene identifiers with gene names. The resulting file will be saved in the `Data_matrix_Smart_seq2` folder (e.g., `matrix_1_Deng.tab`). If matrices already exist, enter `0` when asked for the number of files to process.

2. **Creation of objects (AnnData)**
    - Displays the list of matrices.
    - The script asks for the number of objects to create from the matrices.
    - Creates objects in the `Objects` folder (e.g., `object_1_Deng_ori.h5ad`) and generates a violin plot in the `Plots` folder (e.g., `Plots_1_Deng`) using the `create_anndata_object` function from the `create_anndata_object.py` script.

  ![anndata object](Images/anndata_object_image.png)

### 10X (answer 2 or 3)

1. **Create objects**
    - Displays the list of 10X folders.
    - The script asks for the number of objects to create.
    - Enter the file names one by one.
    - Creates objects using the `create_anndata_object` function (similar to Smart-seq2).

### Filtering Data

1. **Ask if the user wants to filter objects**
    - If yes, displays the list of objects.
    - The script asks for the number of objects to filter.
    - Enter the object names one by one.
    - Executes the `apply_filters` function from the `apply_filters.py` script:
        - Displays the number of cells and genes before filtering.
        - Requests filter parameters (minimum genes per cell, minimum cells per gene, and maximum percentage of mitochondrial genes).
        - Displays the number of cells and genes after filtering and saves a violin plot with the filter parameters in the name.
        - Asks if the user wants to save the filtered object (e.g., `object_1_Deng_filtered.h5ad`).
> [!WARNING]
> If a filtered object already exists, it will be overwritten by the new one.

### Merging Objects (answer 4 or yes at the end of object creation)

1. **Ask if the user wants to merge objects**
    - If Smart-seq2 and/or 10X were processed (answer 1, 2, or 3), ask if the user wants to merge objects. Otherwise, proceed to merge only existing objects (answer 4).
    - Displays the list of objects.
    - The script asks for the number of objects to merge.
    - Enter the object names one by one.
    - Merges objects using the `merge` function from the `merge.py` script and saves the merged object with the names of all included objects in the filename (e.g., `object_1_Deng_filtered_2_Chen_filtered_3_Mohammed_filtered_4_Nowo_merged_6_Posfai_filtered_10_Arg_filtered.h5ad`).
  

## Notes

- Ensure all necessary scripts (`generation_matrix.py`, `create_anndata_object.py`, `apply_filters.py`, `merge.py` and `main.py`) are available in the script directory.
- Review filter parameters and merged object names to avoid overwriting important data.

