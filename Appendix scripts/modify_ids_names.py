import pandas as pd

name_file = ''  # The input file
name_output = ''  # Name for the output file

# Read the matrix data from a tab-separated file into a DataFrame
df = pd.read_table(name_file, sep = '\t', header = [0])

# Read gene names from 'gene_names.tab' into a DataFrame
df_gene_name = pd.read_table('../gene_names.tab', header = None, names = ['GeneID', 'gene_name'])

# Merge the gene names DataFrame with the matrix DataFrame on 'GeneID'
merged_df = pd.merge(df_gene_name, df, on = 'GeneID')

# Delete the 'GeneID' column from the merged DataFrame
del merged_df['GeneID']

# Save the merged DataFrame to a tab-separated file
merged_df.to_csv('Data/' + name_output, sep = '\t', index = False)

print('Done')
