import os, subprocess

dataDirectory = os.path.join(os.path.dirname(os.path.abspath(__file__)),"data")
geneDesignationsToParseFilePath = os.path.join(dataDirectory,"C_elegans_gene_designations.tsv")
geneDesignationsParsedFilePath = os.path.join(dataDirectory,"C_elegans_gene_designations.bed")

with open(geneDesignationsToParseFilePath, 'r') as geneDesignationsToParseFile:
    with open(geneDesignationsParsedFilePath, 'w') as geneDesignationsParsedFile:

        geneDesignationsToParseFile.readline()

        for line in geneDesignationsToParseFile:

            choppedUpLine = line.strip().split('\t')

            chromosome = choppedUpLine[2]
            # I'm not actually sure if the designations in this file are 0- or 1-based.  
            # I'm just assuming they are 1-based for now.,
            startPos0Based = str(int(choppedUpLine[4]) - 1)
            EndPos1Based = choppedUpLine[5]
            strand = choppedUpLine[3]

            geneDesignationsParsedFile.write('\t'.join((chromosome,startPos0Based,EndPos1Based,strand)) + '\n')

subprocess.run(" ".join(("sort","-k1,1","-k2,2n",geneDesignationsParsedFilePath,
                         "-o",geneDesignationsParsedFilePath)), shell = True, check = True)