##### Take genes names and genes ids

# Delete spaces (easier for awk)
sed -i 's/ //g' mm10_all_min10_extended.gtf

# id_gene recovery and gene name in column 9 of gtf, separated by a tab
file=/mnt/c/Users/learogue/Documents/mm10_all_min10_extended.gtf
awk '{print $9}' $file | awk -F ';' '{print $1 "\t" $3}'| sed 's/gene_id//' | sed 's/gene_name//' | tr -d '"' > gene_names.tab

# Delete duplicates
sort -u -o /mnt/c/Users/learogue/Documents/gene_names.tab  /mnt/c/Users/learogue/Documents/gene_names.tab