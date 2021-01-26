import os
from nucperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog


# Combines data from the filtered genes file path and the unfiltered gene designations file path to 
# produce a file of filtered gene designations.
def filterGeneDesignationsByExpression(filteredGenesFilePath: str, unfilteredGeneDesignationsFilePath: str):

    print("Working with:",os.path.basename(filteredGenesFilePath),
          "and",os.path.basename(unfilteredGeneDesignationsFilePath))

    assert filteredGenesFilePath.endswith("_highly_expressed_genes.tsv"), "Unexpected file ending"

    # Get the relevant paths for the input and output files.
    dataDirectory = os.path.dirname(filteredGenesFilePath)
    baseFilteredFileName = os.path.basename(filteredGenesFilePath).rsplit("s.tsv",1)[0]
    filteredGeneDesignationsFilePath = os.path.join(dataDirectory, baseFilteredFileName+"_designations.bed")

    # Get a hash of the highly expressed genes
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


def main():

    # Create a simple dialog for selecting the gene expression file.
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data"))
    dialog.createFileSelector("Filtered Genes File", 0, ("Tab Separated Values File", ".tsv"))
    dialog.createFileSelector("Unfiltered Gene Designations File", 1, ("Bed File", ".bed"))

    dialog.mainloop()

    if dialog.selections is None: quit()

    # Retrieve the file paths to the peak files and pass them to the main funciton.  You know, like, the "main" function.
    # Not this one since it's actually the main funciton, but the one that actually has the main functionality ohmygoodnesswhydoIcodelikethis
    filterGeneDesignationsByExpression(dialog.selections.getIndividualFilePaths()[0], 
                                       dialog.selections.getIndividualFilePaths()[1])


if __name__ == "__main__": main()