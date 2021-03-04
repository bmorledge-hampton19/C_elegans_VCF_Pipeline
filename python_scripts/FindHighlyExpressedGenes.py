import os, subprocess
from typing import List, Dict
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog

# Given a file of RPKM/FPKM values, pull out the most highly expressed genes based on a given percent cutoff value.
# Filtering may be specified to only use data from a specific tissue/cell type.
def findHighlyExpressedGenes(geneExpressionFilePath, tissueFiltering: str, percentCutoff = 25):

    print("Filtering in",os.path.basename(geneExpressionFilePath),"with a",percentCutoff,"percent cutoff.")
    print("Tissue filtering is set to:",tissueFiltering)
    print()

    # Get the relevant paths for the input and output files.
    dataDirectory = os.path.dirname(geneExpressionFilePath)
    geneExpressionFileName = os.path.basename(geneExpressionFilePath).rsplit('_gene_expression.txt',1)[0]
    if tissueFiltering != "Any": geneExpressionFileName += '_' + tissueFiltering.lower().replace(' ', '_')
    averageGeneExpressionFilePath = os.path.join(dataDirectory, geneExpressionFileName + "_average_FPKM.tsv")
    filteredGenesFilePath = os.path.join(dataDirectory, geneExpressionFileName + "_highly_expressed_genes.tsv")

    with open(geneExpressionFilePath, 'r') as geneExpressionFile:
        with open(averageGeneExpressionFilePath, 'w') as averageGeneExpressionFile:

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

            # Populate the dictionary of raw RPKM values if the data set meets the conditions.

            for line in geneExpressionFile:

                choppedUpLine = line.strip().split('\t')

                # Filter as necessary.
                if tissueFiltering == "Sperm Only":
                    if not "SPERM" in choppedUpLine[4].upper(): continue
                elif tissueFiltering == "Oocyte Only":
                    if not "OOCYTE" in choppedUpLine[4].upper(): continue
                elif tissueFiltering == "Germ Line Only":
                    if not "GERM LINE" in choppedUpLine[4].upper(): continue
                elif tissueFiltering == "Any Germ Line Association":
                    if not ("SPERM" in choppedUpLine[4].upper() or "OOCYTE" in choppedUpLine[4].upper() or
                            "GERMLINE" in choppedUpLine[4].upper()): continue
                else: assert tissueFiltering == "Any", "Unrecognized tissue filtering option: " + tissueFiltering

                for i,RPKMValue in enumerate(choppedUpLine[7:]):
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


def main():

    # Create a simple dialog for selecting the gene expression file.
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data"))
    dialog.createFileSelector("Gene Expression File", 0, ("Text File", ".txt"))
    dialog.createDropdown("Cell/Tissue Type Filtering", 1, 0, ("Any","Sperm Only", "Oocyte Only", 
                                                               "Germ Line Only", "Any Germ Line Association"))
    dialog.mainloop()

    if dialog.selections is None: quit()

    # Retrieve the selections and pass the relevant arguments to the primary function.
    findHighlyExpressedGenes(dialog.selections.getIndividualFilePaths()[0], dialog.selections.getDropdownSelections()[0])


if __name__ == "__main__": main()