# TRUNCATING PART
mutDATA <- read.csv (file = 'SOURCE_ccletrimmed_namechanged_HRT18removed_HEC1Bremoved.csv')
colnames(mutDATA)[8] <- "case_id"
TRUNCATIONmutations = subset(mutDATA, mutation_type == 'Nonsense_Mutation') 

# first merge genes motif data with truncation mutations
merged <- merge(TRUNCATIONmutations, genbankkinases, by = 'Gene', nomatch = 0)


# will match ref amino acid as well as codon <endCatalytic) - AS A CHECKING MEASURE LOOK AT HOW THIS CHANGES IF JUST MATCHING BY CODON
# extract Codon / Ref  from amino acid change
merged$Codon <- as.numeric(str_extract(merged$amino_acid_change, "[0-9]+"))
merged$RefAmino <- str_extract(merged$amino_acid_change, "[A-Z]")

# extract predicted amino acid from each protein seq (calculated from Codon)
merged$codonConfirm <- substr(merged$Protein_Seq, merged$Codon, merged$Codon)

#to ensure cases with more than one truncation mutation of the same gene do not skew results remove duplicate of case_id + transcript (not gene as this may remove incorrect transcript) - potentially this could leave an invalid mutation and take out a valid - but will manually check an invalids to make sure there isnt a valid - rank with smallest first to try and mitigate
mergedNODUPS <- merged[order(merged$Codon),]
mergedNODUPS$noMULTIPLEhits <- paste(mergedNODUPS$f, mergedNODUPS$case_id)
mergedNODUPS$transcriptDUP = !duplicated(mergedNODUPS$noMULTIPLEhits)
mergedNODUPS <- subset(mergedNODUPS, transcriptDUP == 'TRUE')

#merge by Codon smaller than endCAT and Ref Amino acid is the same as locAmino
truncSTRINGENT <- subset (mergedNODUPS, (Codon <= endCAT) & (RefAmino == codonConfirm))

#rank matches by ascending endCAT result so that the shortest matching transcript is retained
truncSTRINGENTSHORTEST <- truncSTRINGENT[order(truncSTRINGENT$endCAT),]

#to filter out duplicates of Gene, Codon and Sample (this should leave the shortest transcript)
truncSTRINGENTSHORTEST$duplicateREF <- paste(truncSTRINGENTSHORTEST$Gene, truncSTRINGENTSHORTEST$Codon, truncSTRINGENTSHORTEST$case_id, sep ='_')
truncSTRINGENTSHORTEST$mergedDUP= !duplicated(truncSTRINGENTSHORTEST$duplicateREF)
genefreqSHORTEST <- subset(truncSTRINGENTSHORTEST, mergedDUP == 'TRUE')
genefreqSHORTEST$transcriptFREQ <- paste(genefreqSHORTEST$f, genefreqSHORTEST$endCAT, sep = '_')
freqSTRINGENTSHORTEST <-table(genefreqSHORTEST$transcriptFREQ)
write.csv (freqSTRINGENTSHORTEST, file ='TEMP_freqTRANSCIPTSSHORTESTccleALL.csv')

#to length correct
freqSHORTEST <- read.csv ('TEMP_freqTRANSCIPTSSHORTESTccleALL.csv')
colnames(freqSHORTEST)[2] <- "transcriptFREQ"
#extract gene and endCAT from each freq name
freqSHORTEST$Gene <- gsub    ( "[-].*$"     , "", freqSHORTEST$transcriptFREQ)
freqSHORTEST$endCAT <-  gsub ( ".*[_]"    , "", freqSHORTEST$transcriptFREQ)
freqSHORTEST$endCAT <- as.numeric(freqSHORTEST$endCAT)
#calculate length corrected freq for each transcript entry
freqSHORTEST$oneoverendCAT <- 1 / freqSHORTEST$endCAT
freqSHORTEST$score <- freqSHORTEST$Freq * freqSHORTEST$oneoverendCAT

#add scores for each gene together to get total(using dplyr)
freqADDEDSHORTEST <- ddply(freqSHORTEST,"Gene",numcolwise(sum))

write.csv (freqADDEDSHORTEST, file ='TEMP_lengthcorrectedfreq_ccleALLSHORTEST.csv')

#same for longest transcript
truncSTRINGENTLONGEST <- truncSTRINGENT[order(-truncSTRINGENT$endCAT),]
truncSTRINGENTLONGEST$duplicateREF <- paste(truncSTRINGENTLONGEST$Gene, truncSTRINGENTLONGEST$Codon, truncSTRINGENTLONGEST$case_id, sep ='_')
truncSTRINGENTLONGEST$mergedDUP= !duplicated(truncSTRINGENTLONGEST$duplicateREF)
genefreqLONGEST <- subset(truncSTRINGENTLONGEST, mergedDUP == 'TRUE')
genefreqLONGEST$transcriptFREQ <- paste(genefreqLONGEST$f, genefreqLONGEST$endCAT, sep = '_')
freqSTRINGENTLONGEST <-table(genefreqLONGEST$transcriptFREQ)
write.csv (freqSTRINGENTLONGEST, file ='TEMP_freqTRANSCIPTSLONGESTccleALL.csv')
freqLONGEST <- read.csv ('TEMP_freqTRANSCIPTSLONGESTccleALL.csv')
colnames(freqLONGEST)[2] <- "transcriptFREQ"
freqLONGEST$Gene <- gsub    ( "[-].*$"     , "", freqLONGEST$transcriptFREQ)
freqLONGEST$endCAT <-  gsub ( ".*[_]"    , "", freqLONGEST$transcriptFREQ)
freqLONGEST$endCAT <- as.numeric(freqLONGEST$endCAT)
freqLONGEST$oneoverendCAT <- 1 / freqLONGEST$endCAT
freqLONGEST$score <- freqLONGEST$Freq * freqLONGEST$oneoverendCAT
freqADDEDLONGEST <- ddply(freqLONGEST,"Gene",numcolwise(sum))
write.csv (freqADDEDLONGEST, file ='TEMP_lengthcorrectedfreq_ccleALLLONGEST.csv')

#make table to bind to CCLE genefreqSHORTEST to validate all mutations
validateCCLE <- as.data.frame(genefreqSHORTEST$Gene)
colnames(validateCCLE)[1] <- "Gene"
validateCCLE$transcript <- genefreqSHORTEST$f
validateCCLE$mutation <- genefreqSHORTEST$amino_acid_change
validateCCLE$study <- genefreqSHORTEST$case_id
validateCCLE$case_id <- genefreqSHORTEST$case_id
validateCCLE$codon <- genefreqSHORTEST$Codon
validateCCLE$endCAT <- genefreqSHORTEST$endCAT


#This script doesn't have a final output!!

