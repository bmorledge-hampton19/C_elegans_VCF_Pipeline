import os
from typing import List, Dict
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog


# Given a list of gene designations file paths, check the amount of overlap between genes in each unique pair.
def checkOverlap(geneDesignationsFilePaths, intersectOutputFilePath, unionOutputFilePath):

    print("Working with files:")
    for geneDesignationsFilePath in geneDesignationsFilePaths:
        print('\t' + os.path.basename(geneDesignationsFilePath))
    print()

    # Read all of the given files into a series of dictionaries for easy lookup later.
    geneNameDicts: List[Dict] = list()
    for geneDesignationsFilePath in geneDesignationsFilePaths:

        geneNameDict = dict()

        with open(geneDesignationsFilePath, 'r') as geneDesignationsFile:
            for line in geneDesignationsFile: geneNameDict[line.strip().split('\t')[0]] = None

        geneNameDicts.append(geneNameDict)

        print(len(geneNameDict),"genes found in",os.path.basename(geneDesignationsFilePath))
    print()

    # Check each unique pairing of gene name dictionaries for the amount of overlap present.
    for i in range(len(geneNameDicts)):
        for j in range(i+1, len(geneNameDicts)):

            overlappingGeneCount = 0
            for gene in geneNameDicts[i]:
                if gene in geneNameDicts[j]: overlappingGeneCount += 1
            totalGenes = len(geneNameDicts[i]) + len(geneNameDicts[j]) - overlappingGeneCount

            print("Overlap between",os.path.basename(geneDesignationsFilePaths[i]),
                  "and",os.path.basename(geneDesignationsFilePaths[j]))
            print('\t',overlappingGeneCount,'/',totalGenes, sep = '')
            print()


    # Check total overlap.
    overlappingGenes = dict()
    geneUnion = geneNameDicts[0].copy()
    for gene in geneNameDicts[0]:
        overlappingGenes[gene] = True

    for i in range(1, len(geneNameDicts)):
        geneUnion.update(geneNameDicts[i])
        for gene in overlappingGenes:
            if gene not in geneNameDicts[i]: overlappingGenes[gene] = False

    # Count up the overlapping genes.  
    totalOverlapCount = 0
    for gene in overlappingGenes:
        if overlappingGenes[gene]: 
            totalOverlapCount += 1

    # Output the intersecting genes to the intersection file if requested.
    if intersectOutputFilePath is not None:
        with open(intersectOutputFilePath, 'w') as intersectOutputFile:
            for gene in overlappingGenes:
                if overlappingGenes[gene]: intersectOutputFile.write(gene + '\t' + "NA" + '\n')

    # Output the union genes to the intersection file if requested.
    if unionOutputFilePath is not None:
        with open(unionOutputFilePath, 'w') as unionOutputFile:
            for gene in geneUnion: unionOutputFile.write(gene + '\t' + "NA" + '\n')

    print("Total Overlap: ", totalOverlapCount, '/', len(geneUnion), sep = '')


def main():

    # Create a simple dialog for selecting the gene designation files.
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data"))
    dialog.createMultipleFileSelector("Filtered Genes File", 0, "highly_expressed_genes.tsv", 
                                      ("Tab Separated Values File", ".tsv"))
    intersectDialog = dialog.createDynamicSelector(1, 0)
    intersectDialog.initCheckboxController("Output intersection to file")
    intersectFileSelector = intersectDialog.initDisplay(1, "intersect")
    intersectFileSelector.createFileSelector("Intersection Output", 1, ("Tab Separated Values File", ".tsv"), newFile = True)
    intersectDialog.initDisplayState()

    unionDialog = dialog.createDynamicSelector(2, 0)
    unionDialog.initCheckboxController("Output union to file")
    unionFileSelector = unionDialog.initDisplay(1, "union")
    unionFileSelector.createFileSelector("Union Output", 1, ("Tab Separated Values File", ".tsv"), newFile = True)
    unionDialog.initDisplayState()

    dialog.mainloop()

    if dialog.selections is None: quit()

    # Retrieve the selections and pass the relevant arguments to the primary function.
    intersectOutputFilePath = None
    unionOutputFilePath = None

    if intersectDialog.getControllerVar():
        intersectOutputFilePath = dialog.selections.getIndividualFilePaths("intersect")[0]
    if unionDialog.getControllerVar():
        unionOutputFilePath = dialog.selections.getIndividualFilePaths("union")[0]

    checkOverlap(dialog.selections.getFilePathGroups()[0], intersectOutputFilePath, unionOutputFilePath)


if __name__ == "__main__": main()