import os, subprocess

def main():
    
    # Specify relevant paths to be accessed later.
    dataDirectory = os.path.join(os.path.dirname(os.path.abspath(__file__)),"data")
    SNVDirectory = os.path.join(dataDirectory, "Filtered_VCFS_all_text", "SNV")
    sampleInfoFilePath = os.path.join(dataDirectory, "sample_info.csv")

    genotypes = set() # A hash to keep track of all unique genotypes.
    genotypesBySample = dict() # A dictionary of genotypes with sample IDs as keys.
    mutagens = set() # A hash to keep track of all unique mutagens used.
    mutagensBySample = dict() # A dictionary of mutagens with sample IDs as keys


    # Get information on samples from the csv file.
    with open(sampleInfoFilePath, 'r') as sampleInfoFile:

        # Skip the header lines.
        sampleInfoFile.readline(); sampleInfoFile.readline()

        # Specify which mutagen identifiers are acceptable for the scope of our analysis
        # NO LONGER SPECIFYING HERE.  TAKE ALL THE MUTAGENS!!
        # acceptableMutagens = ("MMS", "EMS", "DMS", "My_Other_Mutagen")

        # Read through every line in the sample info file, pulling out the relevant samples and associating them with their genotype.
        for line in sampleInfoFile:

            choppedUpLine = line.strip().split(',')
            sample = choppedUpLine[0]
            genotype = choppedUpLine[1]
            genotype = genotype.replace(' ','').replace('(','_').replace(')','')
            mutagen = choppedUpLine[4]
            mutagen = mutagen.replace('/','+')
            

            assert sample not in genotypesBySample, "Duplicate sample found: " + sample
            genotypesBySample[sample] = genotype
            genotypes.add(genotype)
            mutagensBySample[sample] = mutagen
            mutagens.add(mutagen)


    mutagenFilePathsByGenotype = dict() # A dictionary of bed file paths with genotypes as keys.
    mutagenFilesByGenotype = dict() # A dictionary of bed files with genotypes as keys.

    # For every mutagen, create a file to store mutations in.
    # Genotype information is stored in each mutation line as the "cohort"
    for mutagen in mutagens:

        mutagenFileDirectory = os.path.join(dataDirectory, "C_elegans_bed_SNVs_for_nucperiod", mutagen)
        if not os.path.exists(mutagenFileDirectory): os.makedirs(mutagenFileDirectory)
        mutagenFilePath = os.path.join(mutagenFileDirectory, mutagen+"_custom_input.bed")
        mutagenFile = open(mutagenFilePath, 'w')
        mutagenFilePathsByGenotype[mutagen] = mutagenFilePath
        mutagenFilesByGenotype[mutagen] = mutagenFile
        

    # Loop through the SNV directory for each of the vcf files within.
    for VCFFileName in os.listdir(SNVDirectory):

        # Make sure this file has an acceptable sample ID.  (If not, something is wrong...)
        sampleName = VCFFileName.split('_')[0]
        assert sampleName in mutagensBySample, "Unknown sample: " + sampleName + " found."
        print("Writing mutations for sample:", sampleName)

        # Read through the VCF file, counting every mutation line.  
        # (This section could be modified in the future to retrieve more specific information on mutations if need be.)
        with open(os.path.join(SNVDirectory,VCFFileName), 'r') as VCFFile:

            mutagen = mutagensBySample[sampleName]
            genotype = genotypesBySample[sampleName]

            for line in VCFFile:

                # Skip header lines.
                if line.startswith('#'): continue

                # Format the mutation for bed file format, and output it to the respective bed file.
                choppedUpLine = line.strip().split('\t')
                chromosome = "chr"+choppedUpLine[0]
                mutationPos0Based = str(int(choppedUpLine[1]) - 1)
                mutationPos1Based = choppedUpLine[1]
                referenceBase = choppedUpLine[3]
                alternateBase = choppedUpLine[4]
                strand = '+'

                mutagenFilesByGenotype[mutagen].write('\t'.join( (chromosome,mutationPos0Based,mutationPos1Based,
                                                                    referenceBase, alternateBase, strand, genotype) ) + '\n' )

    # Close and sort the bed mutation files
    for mutagen in mutagens:
        mutagenFilesByGenotype[mutagen].close()
        subprocess.run(" ".join(("sort","-k1,1","-k2,2n",mutagenFilePathsByGenotype[mutagen],
                                "-o",mutagenFilePathsByGenotype[mutagen])), 
                       shell = True, check = True)

if __name__ == "__main__": main()