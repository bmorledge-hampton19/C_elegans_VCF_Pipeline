# This script takes a mutation file and the coordinates for gene regions (and the transcribed strand) 
# and calculates how many mutations occured transcribed vs. non-transcribed strands.
# NOTE:  Both input files must be sorted for this script to run properly. 
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))

import os, warnings
from typing import List
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (Metadata, generateFilePath, getDataDirectory,
                                                                  DataTypeStr, getContext)

class MutationData:

    def __init__(self, line, acceptableChromosomes):

        # Read in the next line.
        choppedUpLine = line.strip().split()

        # Assign variables
        self.chromosome = choppedUpLine[0] # The chromosome that houses the mutation.
        self.position = int(choppedUpLine[1]) # The position of the mutation in its chromosome. (0 base)
        self.context = choppedUpLine[3] + '>' + choppedUpLine[4] # The context of the given mutation.
        self.strand = choppedUpLine[5] # Either '+' or '-' depending on which strand houses the mutation.
        self.strandMatchesTS = None # Whether or not the strand of the mutation matches the transcribed strand in
                                    # the gene region encompassing the mutation (NoneType if intergenic or ambiguous)

        # Make sure the mutation is in a valid chromosome.
        if not self.chromosome in acceptableChromosomes:
            raise ValueError(choppedUpLine[0] + " is not a valid chromosome for the mutation trinuc file.")


# Contains data on a single gene position obtained by reading the next available line in a given file.
class GeneData:

    def __init__(self, line):

        # Read in the next line.
        choppedUpLine = line.strip().split()

        self.chromosome = choppedUpLine[0]
        self.startPos = int(choppedUpLine[1]) # 0 base
        self.endPos = int(choppedUpLine[2]) - 1 # Still 0 base 
        reverser = {'+':'-', '-':'+'}
        self.transcribedStrand = reverser[choppedUpLine[5]] # The transcribed strand (reversed because the coding strand is given)


# Uses the given gene positions file and mutation file to count the number of mutations 
# in transcribed, non-transcribed and intergenic regions.  
# Generates a new file to store these results.
# NOTE:  It is VITAL that both files are sorted, first by chromosome number and then by starting coordinate.
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))
#        This code is pretty slick, but it will crash and burn and give you a heap of garbage as output if the inputs aren't sorted.
class CountsFileGenerator:

    def __init__(self, mutationFilePath, genePositionsFilePath, 
                 transcribedRegionMutationCountsFilePath, acceptableChromosomes):

        # Open the mutation and gene positions files to compare against one another.
        self.mutationFile = open(mutationFilePath, 'r')
        self.genePosFile = open(genePositionsFilePath,'r')

        # Store the other arguments passed to the constructor
        self.acceptableChromosomes = acceptableChromosomes
        self.transcribedRegionMutationCountsFilePath = transcribedRegionMutationCountsFilePath

        # Dictionaries holding the number of mutations found for each context.
        # Key is in format "NNN>N"
        self.transcribedRegionMutationCounts = dict()
        self.nontranscribedRegionMutationCounts = dict() 
        self.intergenicAndAmbiguousMutationCounts = dict()

        # Keeps track of mutations that matched to a gene to check for overlap.
        self.mutationsInPotentialOverlap: List[MutationData] = list()

        # The mutation and gene currently being investigated.
        self.currentMutation: MutationData = None
        self.currentGene: GeneData = None


    # Reads in the next mutation from the mutation data into currentMutation
    def readNextMutation(self) -> MutationData:

        # Read in the next line.
        nextLine = self.mutationFile.readline()

        # Check if EOF has been reached.
        if len(nextLine) == 0: 
            self.currentMutation = None
        # Otherwise, read in the next mutation.
        else:
            self.currentMutation = MutationData(nextLine, self.acceptableChromosomes)
            # Assign each mutation to intergenic/ambiguous by default
            self.intergenicAndAmbiguousMutationCounts[self.currentMutation.context] = \
                self.intergenicAndAmbiguousMutationCounts.setdefault(self.currentMutation.context,0) + 1

    
    # Reads in the next gene from the genePos file into currentGene
    def readNextGene(self) -> GeneData:

        # Read in the next line.
        nextLine = self.genePosFile.readline()

        # Check if EOF has been reached.
        if len(nextLine) == 0: 
            self.currentGene = None
        # Otherwise, read in the next mutation.
        else:
            self.currentGene = GeneData(nextLine)

        # Check for mutations in overlapping regions between this gene and the last one.
        if self.currentGene is not None: self.checkMutationsInOverlap() 


    # Takes a mutation object and gene object which have unequal chromosomes and read through data until they are equal.
    def reconcileChromosomes(self):
        
        chromosomeChanged = False # A simple flag to determine when to inform the user that a new chromosome is being accessed

        # Until the chromosomes are the same for both mutations and genes, read through the one with the eariler chromosome.
        while (self.currentMutation is not None and self.currentGene is not None and 
               self.currentMutation.chromosome != self.currentGene.chromosome):
            chromosomeChanged = True
            if self.currentMutation.chromosome < self.currentGene.chromosome: self.readNextMutation()
            else: self.readNextGene()

        if chromosomeChanged and self.currentGene is not None and self.currentMutation is not None: 
            print("Counting in",self.currentGene.chromosome)


    # Determines whether or not the current mutation is past the range of the current gene.
    def isMutationPastGene(self):

        if self.currentMutation is None:
            return True
        elif self.currentMutation.position > self.currentGene.endPos:
            return True
        elif not self.currentMutation.chromosome == self.currentGene.chromosome:
            return True
        else: 
            return False


    # A function which checks to see if the current mutation falls within the current gene and if it does, adds it to the list.
    # (Assumes that the current mutation is not beyond the gene because this is checked first by isMutationPastGene.)
    def addMutationIfInGene(self):

        if self.currentMutation.position >= self.currentGene.startPos:

            # Move this mutation from intergenic counts to the TS/NTS counts
            self.intergenicAndAmbiguousMutationCounts[self.currentMutation.context] -= 1
            if self.currentMutation.strand == self.currentGene.transcribedStrand: 
                self.transcribedRegionMutationCounts[self.currentMutation.context] = \
                    self.transcribedRegionMutationCounts.setdefault(self.currentMutation.context, 0) + 1
                self.currentMutation.strandMatchesTS = True
            else: 
                self.nontranscribedRegionMutationCounts[self.currentMutation.context] = \
                    self.nontranscribedRegionMutationCounts.setdefault(self.currentMutation.context, 0) + 1
                self.currentMutation.strandMatchesTS = False
            
            # Add the mutation to the list of mutations in the current nucleosome
            self.mutationsInPotentialOverlap.append(self.currentMutation)


    # Check to see if any previous mutations called for previous genes are present in the current gene due to overlap.
    def checkMutationsInOverlap(self):    

        # First, get rid of any mutations that fall before the start position of the new gene.
        self.mutationsInPotentialOverlap = [mutation for mutation in self.mutationsInPotentialOverlap 
                                            if mutation.position >= self.currentGene.startPos and 
                                            mutation.chromosome == self.currentGene.chromosome]

        # Next, check all remaining mutations to see if their previous TS/NTS assignment matches with the new gene.
        for mutation in self.mutationsInPotentialOverlap:
            
            # Is the mutation within the upper bound of the gene (endPos)?
            # Does the call for the previous gene(s) match this gene?
            if (mutation.position <= self.currentGene.endPos and 
                mutation.strandMatchesTS != (mutation.strand == self.currentGene.transcribedStrand)):

                # If not, switch the mutation to count for ambiguous.
                if mutation.strandMatchesTS: self.transcribedRegionMutationCounts[mutation.context] -= 1
                else: self.nontranscribedRegionMutationCounts[mutation.context] -= 1
                self.intergenicAndAmbiguousMutationCounts[mutation.context] += 1
                mutation.strandMatchesTS = None

        # Remove any mutations that were found to be ambiguous.
        self.mutationsInPotentialOverlap = [mutation for mutation in self.mutationsInPotentialOverlap 
                                            if mutation.strandMatchesTS is not None]


    # Assign all mutations to either the TS, NTS, or intergenic/ambiguous bins based on the genePos file.
    # (Further bins results by mutation context.)
    def count(self):
        # Get data on the first mutation and gene and reconcile their chromosomes if necessary to start things off.
        # If either the mutation file or gene file is empty, make sure to bypass the check.
        self.readNextMutation()
        self.readNextGene()
        if self.currentMutation is None or self.currentGene is None:
            warnings.warn("Empty Mutation or Gene Positions file.  Output will most likely be unhelpful.")
        elif self.currentGene.chromosome == self.currentMutation.chromosome: 
            print("Counting in",self.currentGene.chromosome)
        else: self.reconcileChromosomes()

        # The core loop goes through each gene one at a time and checks mutation positions against it until 
        # one exceeds its rightmost position or is on a different chromosome (or mutations are exhausted).  
        # Then, the next gene is checked, then the next, etc. until no genes remain.
        while self.currentGene is not None:

            # Read mutations until the mutation is past the range of the current gene.
            while not self.isMutationPastGene():

                # Check and see if we need to add the mutation to our lists.
                self.addMutationIfInGene()
                #Get data on the next mutation.
                self.readNextMutation()

            # Read in a new gene.
            self.readNextGene()

            # Reconcile the mutation data and gene data to be sure
            # that they are looking at the same chromosome for the next iteration
            self.reconcileChromosomes()


    def writeResults(self):
        # Write the results to the output file.
        with open(self.transcribedRegionMutationCountsFilePath,'w') as TSMutationCountsFile:
            
            # Write the headers to the file.
            TSMutationCountsFile.write('\t'.join(("Mutation_Context", "Mutant_Base", "TS_Counts",
                                                  "NTS_Counts","Intergenic_And_Ambiguous_Counts",
                                                  "NTS_to_TS_Ratio")) + '\n')

            # Write the counts.
            for context in sorted(self.intergenicAndAmbiguousMutationCounts.keys()):

                TSCounts = self.transcribedRegionMutationCounts.setdefault(context,0)
                NTSCounts = self.nontranscribedRegionMutationCounts.setdefault(context,0)
                IACounts = self.intergenicAndAmbiguousMutationCounts[context]

                # Don't divide by 0 >:(
                if TSCounts != 0: NTSOverTS = NTSCounts/TSCounts
                else: NTSOverTS = "NA"

                TSMutationCountsFile.write('\t'.join((context.split('>')[0],context.split('>')[1],
                                                      str(TSCounts), str(NTSCounts), str(IACounts), str(NTSOverTS))) + '\n')
            

# Main functionality starts here.
def countInTranscribedRegions(mutationFilePaths, genePositionsFilePath):

    transcribedRegionMutationCountsFilePaths = list() # A list of paths to the output files generated by the function

    # Loop through each given mutation file path, creating a corresponding transcribed region mutation count file for each.
    for mutationFilePath in mutationFilePaths:

        print("\nWorking with",os.path.split(mutationFilePath)[1])

        # Make sure we have the expected file type.
        if not DataTypeStr.mutations in os.path.basename(mutationFilePath): 
            raise ValueError("Mutation file should have \"" + DataTypeStr.mutations + "\" in the name.")
        
        assert getContext(mutationFilePath, True) == 3, ("Expected trinucleotide context mutation file." +
                                                         "Code needs to be modified to accept other formats.")
            
        # Get the list of acceptable chromosomes
        acceptableChromosomes = ("chrI","chrII","chrIII","chrIV","chrV","chrX")
        
        # Get metadata and use it to generate a path to the nucleosome positions file.
        metadata = Metadata(mutationFilePath)

        # Generate the output file path
        transcribedRegionMutationCountsFilePath = generateFilePath(directory = metadata.directory,
                                                                   dataGroup = metadata.dataGroupName, 
                                                                   fileExtension = ".tsv", dataType = "transcriptional_asymmetry")

        # Ready, set, go!
        counter = CountsFileGenerator(mutationFilePath, genePositionsFilePath, 
                                      transcribedRegionMutationCountsFilePath, acceptableChromosomes)
        counter.count()
        counter.writeResults()
        transcribedRegionMutationCountsFilePaths.append(transcribedRegionMutationCountsFilePath)

    return transcribedRegionMutationCountsFilePaths


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Mutation Files:",0,DataTypeStr.mutations + ".bed",("Bed Files",".bed"))
    dialog.createFileSelector("Gene Position File:", 1, ("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    mutationFilePaths = selections.getFilePathGroups()[0] # A list of mutation file paths
    genePositionsFilePath = selections.getIndividualFilePaths()[0] # The gene positions file path

    countInTranscribedRegions(mutationFilePaths, genePositionsFilePath)

if __name__ == "__main__": main()