##########################
# Script 3: Merging the truncating mutation data from the TCGA and CCLE truncating scripts (scripts 1 and 2)
#
# Automatic merge function to combine the TCGA and CCLE trucation mutation data outputted from scripts 1 and 2.
##########################

# Intro text.
cat("A. Hudson, N. Stephenson, C. Li, E. Trotter, A. Fletcher, G. Katrona, P. Bieniasz-Krzywiec, M. Howell, C. Writh, C. Miller, J. Brognard. \n\n",
    "Functional Screening of Cancer Datasets Identifies Mutational Hotspots and Novel Targets in the Tumour-Suppressing Kinome. \n\n", 
    "Script requires stringr and plyr packages, additionally set the working directory to the folder containing the SOURCE.csv files \n\n", 
    sep="")

# merge two SHORTEST truncating scores (using dplyr)
ccletruncSHORTEST <- read.csv (file = "TEMP_lengthcorrectedfreq_ccleALLSHORTEST.csv")
tcgatruncSHORTEST <- read.csv (file = "TEMP_lengthcorrectedfreq_tcgaALLSHORTEST.csv")
combinedTRUNCSHORTEST <- rbind(ccletruncSHORTEST, tcgatruncSHORTEST)
combinedTRUNCscoresSHORTEST <- ddply(combinedTRUNCSHORTEST,"Gene",numcolwise(sum))
write.csv (combinedTRUNCscoresSHORTEST, file ='outputs/mergedTRUNCscoresSHORTEST.csv')

# merge two LONGEST truncating scores (using dplyr)
ccletruncLONGEST <- read.csv (file = 'TEMP_lengthcorrectedfreq_ccleALLLONGEST.csv')
tcgatruncLONGEST <- read.csv (file = 'TEMP_lengthcorrectedfreq_tcgaALLLONGEST.csv')
combinedTRUNCLONGEST <- rbind(ccletruncLONGEST, tcgatruncLONGEST)
combinedTRUNCscoresLONGEST <- ddply(combinedTRUNCLONGEST,"Gene",numcolwise(sum))
write.csv (combinedTRUNCscoresLONGEST, file ='outputs/mergedTRUNCscoresLONGEST.csv')

#get mean of LONGEST and SHORTEST scores
combinedTRUNCscoresMEAN <- combinedTRUNCscoresSHORTEST
combinedTRUNCscoresMEAN$longestSCORE <- as.numeric(combinedTRUNCscoresLONGEST$score)
colnames(combinedTRUNCscoresMEAN)[7] <- "shortestSCORE"
combinedTRUNCscoresMEAN$shortestSCORE <- as.numeric (combinedTRUNCscoresMEAN$shortestSCORE)
combinedTRUNCscoresMEAN$meanSCORE = (combinedTRUNCscoresMEAN$shortestSCORE + combinedTRUNCscoresMEAN$longestSCORE) / 2
write.csv (combinedTRUNCscoresMEAN, file ='outputs/mergedTRUNCscoresMEAN.csv')
SOURCE_Top30 <- combinedTRUNCscoresMEAN[with(combinedTRUNCscoresMEAN, order(-meanSCORE)),]
SOURCE_Top30 <- head(SOURCE_Top30, n=30L)
SOURCE_Top30 <- as.data.frame(SOURCE_Top30$Gene)
colnames(SOURCE_Top30)[1] <- "Gene"
write.csv (SOURCE_Top30, file="SOURCE_Top30.csv")


#make table to bind to CCLE genefreqSHORTEST to validate all mutations
validateTCGA <- as.data.frame(genefreqSHORTEST$Gene)
colnames(validateTCGA)[1] <- "Gene"
validateTCGA$transcript <- genefreqSHORTEST$f
validateTCGA$mutation <- genefreqSHORTEST$amino_acid_change
validateTCGA$study <- genefreqSHORTEST$genetic_profile_id
validateTCGA$case_id <- genefreqSHORTEST$case_id
validateTCGA$codon <- genefreqSHORTEST$Codon
validateTCGA$endCAT <- genefreqSHORTEST$endCAT

#bind to validateCCLE
validateCCLE$study <- NULL
top30check <- rbind (validateTCGA, validateCCLE) 


#check the top30 to make sure the all transcript analysis doesnt give any false positives
top30 <- read.csv ('SOURCE_Top30.csv')
top30check <- merge (top30, top30check, by = 'Gene', nomatch = 0)
write.csv (top30check, file ='top30check.csv')

#these were checked and are all valid - not extract catalytic frag of top 30 - note there are mutliple diff transcripts - use longest one that matches
# use genefreqLONGEST which is all tcga matches ranked by descending size of transcript
#first get rid of duplicates for gene to leave longest transcript
FRAGEXTRACT <- genefreqLONGEST
FRAGEXTRACT$DUP= !duplicated(FRAGEXTRACT$Gene)
FRAGEXTRACT <- subset(FRAGEXTRACT, DUP == 'TRUE')
FRAGEXTRACTtop30 <- merge (top30, FRAGEXTRACT, by = 'Gene', nomatch = 0)
#AMHR2 was changed manually as it is incorrect
write.csv (FRAGEXTRACTtop30, file ='TEMP_check.csv')


