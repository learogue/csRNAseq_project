import glob
import pandas as pd

def generation_matrix(path):
    """
    Create a matrix for files containing one cell in a directory from the same dataset (in Data).

    Input:
        path (str): The counts directory
            Example of a file:
                ENSMUSG00000096126 0
                ENSMUSG00000102735 0
                ENSMUSG00000103922 4
                ENSMUSG00000025903 115

    Output:
        Matrix with all cells of the dataset in Data_matrix_Smart_seq2, format tabular
            gene_name  SRR805347  SRR805348  SRR805349  SRR805350  SRR805351  SRR805352
            Gnai3      1149       983        2556       1171       2123       1119
            Pbsn       0          0          0          0          0          0
    """
    # Create the output file name by replacing '_counts' with '.tab' in the input path
    name_file = 'matrix_' + path.replace('_counts', '.tab')

    # Initializing the empty DataFrame to contain the merged data
    merged_df = pd.DataFrame()

    # Browse each file in the directory and merge data
    for file in glob.glob('Data/' + path + '/*.tabular'):
        # Extract cell name from file path (assuming it follows the format 'Data/path/cell_name.tabular')
        cell_name = file.replace('\\', '/').split('/')[2].replace('.tabular', '')
        
        # Read the count file into a DataFrame
        df = pd.read_table(file, sep = '\t', header = None, names = ['gene_id', cell_name])
        
        # If the merged DataFrame is empty, initialize it with the current DataFrame
        if merged_df.empty:
            merged_df = df
        else:
            # Otherwise, merge the current DataFrame with the merged DataFrame on 'gene_id'
            merged_df = pd.merge(merged_df, df, on = 'gene_id', how = 'outer')

    # Read the file containing gene IDs and names
    df_gene_name = pd.read_table('gene_names.tab', header = None, names = ['gene_id', 'gene_name'])

    # Merge the gene names DataFrame with the merged count DataFrame on 'gene_id'
    merged_df = pd.merge(df_gene_name, merged_df, on = 'gene_id', how = 'right')

    # Delete the 'gene_id' column as it's no longer needed
    del merged_df['gene_id']

    # Save the final merged DataFrame to a tab-separated file in the Data_matrix_Smart_seq2 directory
    merged_df.to_csv('Data_matrix_Smart_seq2/' + name_file, sep = '\t', index = False)

    print('Done')
