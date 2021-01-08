import os, subprocess
from typing import List, Dict

# Get the relevant paths for the input and output files.
dataDirectory = os.path.join(os.path.dirname(os.path.abspath(__file__)),"data")
geneExpressionFilePath = os.path.join(dataDirectory, "C_elegans_gene_expression.txt")
averageGeneExpressionFilePath = os.path.join(dataDirectory, "C_elegans_expressed_genes_average_RPKM.tsv")
filteredGenesFilePath = os.path.join(dataDirectory, "C_elegans_highly_expressed_genes.tsv")

with open(geneExpressionFilePath, 'r') as geneExpressionFile:
    with open(averageGeneExpressionFilePath, 'w') as averageGeneExpressionFile:

        percentCutoff = 25 # What percent of genes to take off the top and consider "highly expressed"
        genes = list() # An ordered list of the genes 
        rawRPKM: Dict[str,List] = dict() # A dictionary containing the list of RPKM values for each gene.

        # Skip the header.
        geneExpressionFile.readline()

        # Populate the ordered list of gene names and create lists in rawRPKM.
        for gene in geneExpressionFile.readline().strip().split('\t')[7:]:
            gene = gene.split('(',1)[1].rsplit(')',1)[0]
            gene = gene.replace(", ", "__")
            genes.append(gene)
            rawRPKM[gene] = list()

        # Populate the dictionary of raw RPKM values
        for line in geneExpressionFile:
            for i,RPKMValue in enumerate(line.strip().split('\t')[7:]):
                #print("Gene:",genes[i])
                if RPKMValue == '' or RPKMValue == "N.A.": continue
                rawRPKM[genes[i]].append(float(RPKMValue))

        # Average the values.  Write the averaged values and their respective genes to the output file.
        for gene in genes:
            if len(rawRPKM[gene]) == 0: averageRPKM = '0'
            else: averageRPKM = str(sum(rawRPKM[gene])/len(rawRPKM[gene]))
            averageGeneExpressionFile.write(gene + '\t' + averageRPKM + '\n')


# Sort the output file in descending order
subprocess.run(" ".join(("sort","-k2,2gr",averageGeneExpressionFilePath,
                         "-o",averageGeneExpressionFilePath)), shell = True, check = True)

# Output the highly expressed genes to a separate file.
geneCountCutoff = str(int(len(genes)*percentCutoff/100))
subprocess.run(" ".join(("head", '-'+geneCountCutoff, averageGeneExpressionFilePath, 
                         '>', filteredGenesFilePath)), shell = True, check = True)