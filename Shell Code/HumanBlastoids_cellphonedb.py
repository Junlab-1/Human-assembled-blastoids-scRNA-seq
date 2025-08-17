from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
result = cpdb_statistical_analysis_method.call(
    cpdb_file_path = '~/software/CellphoneDB/cellphonedb.zip',            # mandatory: CellPhoneDB database zip file.
    meta_file_path = "./HumanBlastoids_meta.txt",              # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = './HumanBlastoids_counts.txt',# mandatory: normalized count matrix.
    counts_data = 'ensembl',   #--counts-data: [ensembl | gene_name | hgnc_symbol] Type of gene identifiers in the counts data
	output_path = './',
    )
