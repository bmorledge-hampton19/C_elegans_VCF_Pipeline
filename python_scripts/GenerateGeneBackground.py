from nucperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from nucperiodpy.helper_scripts.UsefulBioinformaticsFunctions import bedToFasta, reverseCompliment, FastaFileIterator
from typing import List
import os

# Given a genome fasta file path, generates background counts for every trinucleotide sequence for each of the given 
# sets of gene designations.  TS and NTS trinucleotides are counted separately, and ambiguous gene regions (any 
# mixing of '+' and '-' regions) are not counted.
def generateGeneBackground(geneDesignationsFilePaths: List[str], genomeFilePath):

    for geneDesignationsFilePath in geneDesignationsFilePaths:

        print("Generating background for",geneDesignationsFilePath)

        assert geneDesignationsFilePath.endswith("gene_designations.bed"), ("Unexpected file path.  Expected file ending in " +
                                                                            "gene_designations.bed")

        print("Parsing clear ranges...")

        # First, condense all overlapping gene regions and remove any ambiguous regions.
        clearGeneRangesFilePath = geneDesignationsFilePath.rsplit("gene_designations.bed", 1)[0] + "clear_gene_ranges.bed"

        currentGeneRangeChromosome = None
        currentGeneRangeStart = None
        currentGeneRangeEnd = None
        currentGeneRangeStrand = None

        with open(geneDesignationsFilePath, 'r') as geneDesignationsFile:
            with open(clearGeneRangesFilePath, 'w') as clearGeneRangesFile:
                for line in geneDesignationsFile:

                    # Parse out the gene range info from the current line.
                    choppedUpLine = line.strip().split('\t')
                    lineChromosome = choppedUpLine[0]
                    lineGeneStart = int(choppedUpLine[1])
                    lineGeneEnd = int(choppedUpLine[2])
                    lineStrand = choppedUpLine[3]

                    # Unless we are starting a new gene range, check to see if the gene region on this line overlaps with the current one.
                    if currentGeneRangeChromosome is not None:

                        # If they overlap, expand the current range and check to make sure the strands match.
                        if currentGeneRangeChromosome == lineChromosome and lineGeneStart < currentGeneRangeEnd:
                            
                            currentGeneRangeEnd = lineGeneEnd
                            if currentGeneRangeStrand is not None and currentGeneRangeStrand != lineStrand: currentGeneRangeStrand = None

                        # If the don't overlap, check to make sure the strand designation for this region is unambiguous, then write it.
                        # Also, keep in mind to expand the ranges by one bp on either side for trinucleotide context at the borders.
                        else:

                            if currentGeneRangeStrand is not None:
                                clearGeneRangesFile.write('\t'.join((currentGeneRangeChromosome, str(currentGeneRangeStart - 1), 
                                                                     str(currentGeneRangeEnd + 1), currentGeneRangeStrand)) + '\n')
                            
                            # Don't forget to reset the chromosome variable to flag the rest for reassignment!
                            currentGeneRangeChromosome = None


                    # If we are starting to look at a new gene range, assign all the values from this line.
                    if currentGeneRangeChromosome is None:
                        currentGeneRangeChromosome = lineChromosome
                        currentGeneRangeStart = lineGeneStart
                        currentGeneRangeEnd = lineGeneEnd
                        currentGeneRangeStrand = lineStrand

                # Do one last check so we don't miss the last gene range.
                if currentGeneRangeStrand is not None:
                    clearGeneRangesFile.write('\t'.join((currentGeneRangeChromosome, str(currentGeneRangeStart - 1), 
                                                         str(currentGeneRangeEnd + 1), currentGeneRangeStrand)) + '\n')


        # Obtain trinucleotide counts for all the gene ranges.

        # Generate the fasta file...
        print("Generating fasta file...")
        clearGeneRangeFastaFilePath = clearGeneRangesFilePath.rsplit('.',1)[0] + ".fa"
        bedToFasta(clearGeneRangesFilePath, genomeFilePath, clearGeneRangeFastaFilePath)

        # Iterate through the fasta file, counting trinucleotides.
        print("Counting and writing trinucleotides...")
        trinucleotideCountsNTS = dict()
        trinucleotideCountsTS = dict()
        with open(clearGeneRangeFastaFilePath, 'r') as clearGeneRangesFastaFile:
            for fastaEntry in FastaFileIterator(clearGeneRangesFastaFile):

                for i in range(0, len(fastaEntry.sequence) - 2):

                    trinucleotideNTS = fastaEntry.sequence[i:i+3]
                    trinucleotideTS = reverseCompliment(trinucleotideNTS)
                    trinucleotideCountsNTS[trinucleotideNTS] = trinucleotideCountsNTS.setdefault(trinucleotideNTS, 0) + 1
                    trinucleotideCountsTS[trinucleotideTS] = trinucleotideCountsTS.setdefault(trinucleotideTS, 0) + 1

        # Write the background trinucleotide counts to a separate file.
        trinucleotideBackgroundCountsFilePath = clearGeneRangesFilePath.rsplit("clear_gene_ranges.bed",1)[0] + "background_gene_trinuc_counts.tsv"
        
        with open(trinucleotideBackgroundCountsFilePath, 'w') as trinucleotideBackgroundCountsFile:

            # Write the header
            trinucleotideBackgroundCountsFile.write('\t'.join(("Trinucleotide", "NTS_Counts", "TS_Counts")) + '\n')

            # Write the counts.

            for trinucleotide in sorted(trinucleotideCountsNTS.keys() | trinucleotideCountsTS.keys()):

                trinucleotideBackgroundCountsFile.write('\t'.join((trinucleotide, str(trinucleotideCountsNTS.setdefault(trinucleotide, 0)),
                                                                   str(trinucleotideCountsTS.setdefault(trinucleotide, 0)))) + '\n')

        print()


def main():

    # Create a simple dialog for selecting the gene expression file.
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data"))
    dialog.createMultipleFileSelector("Gene Designations Files", 0, "gene_designations.bed", ("Text File", ".bed"))
    dialog.createFileSelector("Genome Fasta File", 1, ("Fasta File", ".fa"))
    dialog.mainloop()

    if dialog.selections is None: quit()

    # Retrieve the selections and pass the relevant arguments to the primary function.
    generateGeneBackground(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0])


if __name__ == "__main__": main()