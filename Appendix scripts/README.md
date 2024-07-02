# Appendix scripts

The script `extraction_gene_name.sh` take gene ids and gene names from gtf format (example below), and generates one tab file with in column 1 the gene ids and column 2 gene names.
```
##gff-version2
##source-versionrtracklayer1.58.0
##date2023-01-16
chr1	havana	exon	3073253	3074322	.	+	.	gene_id"ENSMUSG00000102693";gene_version"1";gene_name"4933401J01Rik";gene_source"havana";gene_biotype"TEC";havana_gene"OTTMUSG00000049935";havana_gene_version"1";transcript_id"ENSMUST00000193812";transcript_version"1";transcript_name"4933401J01Rik-201";transcript_source"havana";transcript_biotype"TEC";havana_transcript"OTTMUST00000127109";havana_transcript_version"1";tag"basic";transcript_support_level"NA";exon_number"1";exon_id"ENSMUSE00001343744";exon_version"1";
chr1	ensembl	exon	3102016	3102125	.	+	.	gene_id"ENSMUSG00000064842";gene_version"1";gene_name"Gm26206";gene_source"ensembl";gene_biotype"snRNA";transcript_id"ENSMUST00000082908";transcript_version"1";transcript_name"Gm26206-201";transcript_source"ensembl";transcript_biotype"snRNA";tag"basic";transcript_support_level"NA";exon_number"1";exon_id"ENSMUSE00000522066";exon_version"1";
```

The script `modify_ids_names.py` modifies gene ids in Smart-seq2 matrix in gene names.
