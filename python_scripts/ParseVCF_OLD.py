import os, subprocess

def main():
    
    # Specify relevant paths to be accessed later.
    projectDirectory = os.path.dirname(os.path.abspath(__file__))
    SNVDirectory = os.path.join(projectDirectory, "Filtered_VCFS_all_text", "SNV")
    sampleInfoFilePath = os.path.join(projectDirectory, "sample_info.csv")

    genotypes = set() # A hash to keep track of all unique genotypes.
    genotypesBySample = dict() # A dictionary of genotypes with sample IDs as keys.

    # Get information on samples from the csv file.
    with open(sampleInfoFilePath, 'r') as sampleInfoFile:

        # Skip the header lines.
        sampleInfoFile.readline(); sampleInfoFile.readline()

        # Specify which mutagen identifiers are acceptable for the scope of our analysis
        acceptableMutagens = ("MMS", "EMS", "DMS", "My_Other_Mutagen")

        # Read through every line in the sample info file, pulling out the relevant samples and associating them with their genotype.
        for line in sampleInfoFile:

            choppedUpLine = line.strip().split(',')
            sample = choppedUpLine[0]
            genotype = choppedUpLine[1]
            mutagen = choppedUpLine[4]

            if mutagen in acceptableMutagens:
                assert sample not in genotypesBySample, "Duplicate sample found: " + sample
                genotypesBySample[sample] = genotype
                genotypes.add(genotype)


    mutationCountsByGenotype = dict() # A dictionary of mutation counts with genotypes as keys.
    genotypeFilePathsByGenotype = dict() # A dictionary of bed file paths with genotypes as keys.
    genotypeFilesByGenotype = dict() # A dictionary of bed files with genotypes as keys.
    
    # Create and open the bed files.
    for genotype in genotypes:
        cElegansSNVDirectory = os.path.join(projectDirectory, "C_elegans_bed_SNVs")
        if not os.path.exists(cElegansSNVDirectory): os.makedirs(cElegansSNVDirectory)
        genotypeFilePath = os.path.join(cElegansSNVDirectory, genotype+".bed")
        genotypeFile = open(genotypeFilePath, 'w')
        genotypeFilePathsByGenotype[genotype] = genotypeFilePath
        genotypeFilesByGenotype[genotype] = genotypeFile
        

    # Loop through the SNV directory for each of the vcf files within.
    for VCFFileName in os.listdir(SNVDirectory):

        # Make sure this file has an acceptable sample ID.  (If not, skip it.)
        sampleName = VCFFileName.split('_')[0]
        if sampleName not in genotypesBySample: continue

        # Read through the VCF file, counting every mutation line.  
        # (This section could be modified in the future to retrieve more specific information on mutations if need be.)
        with open(os.path.join(SNVDirectory,VCFFileName), 'r') as VCFFile:
            mutationCounts = 0

            genotype = genotypesBySample[sampleName]

            for line in VCFFile:

                # Skip header lines.
                if line.startswith('#'): continue

                # Format the mutation for bed file format, and output it to the respective bed file.
                choppedUpLine = line.strip().split('\t')
                chromosome = choppedUpLine[0]
                mutationPos0Based = str(int(choppedUpLine[1]) - 1)
                mutationPos1Based = choppedUpLine[1]
                referenceBase = choppedUpLine[3]
                alternateBase = choppedUpLine[4]

                genotypeFilesByGenotype[genotype].write('\t'.join( (chromosome,mutationPos0Based,mutationPos1Based,
                                                                    referenceBase, alternateBase) ) + '\n' )

                # Count the mutation.
                mutationCounts += 1

            mutationCountsByGenotype[genotype] = mutationCountsByGenotype.setdefault(genotype, 0) + mutationCounts

    totalMutations = 0

    # Close and sort the bed mutation files
    for genotype in genotypes:
        genotypeFilesByGenotype[genotype].close()
        subprocess.run(" ".join(("sort","-k1,1","-k2,2n",genotypeFilePathsByGenotype[genotype],
                                 "-o",genotypeFilePathsByGenotype[genotype])), 
                       shell = True, check = True)


    for genotype in mutationCountsByGenotype:
        print(genotype,": ",mutationCountsByGenotype[genotype], sep = '')

        totalMutations += mutationCountsByGenotype[genotype]

    print("Total Mutations:",totalMutations)

if __name__ == "__main__": main()