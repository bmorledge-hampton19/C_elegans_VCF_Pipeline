import os

# Get the relevant paths for the input and output files.
dataDirectory = os.path.join(os.path.dirname(os.path.abspath(__file__)),"data")
filteredGenesFilePath = os.path.join(dataDirectory, "C_elegans_highly_expressed_genes.tsv")
unfilteredGeneDesignationsFilePath = os.path.join(dataDirectory,"C_elegans_gene_designations.bed")
filteredGeneDesignationsFilePath = os.path.join(dataDirectory, "C_elegans_high_exp_gene_designations.bed")

# Get a list of the highly expressed genes
filteredGenes = dict()
with open(filteredGenesFilePath, 'r') as filteredGenesFile:
    
    for line in filteredGenesFile:
        # NOTE: The double underscore string is used to separate the multiple acceptable names for a single gene.
        for gene in line.strip().split('\t')[0].split('__'): filteredGenes[gene] = None

# Filter the original gene designations file to only include data from the filtered genes.
with open(unfilteredGeneDesignationsFilePath, 'r') as unfilteredGeneDesignationsFile:
    with open(filteredGeneDesignationsFilePath, 'w') as filteredGeneDesignationsFile:

        for line in unfilteredGeneDesignationsFile:
            if line.strip().split('\t')[4] in filteredGenes:
                filteredGeneDesignationsFile.write(line)
