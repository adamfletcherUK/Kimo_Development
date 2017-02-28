
# TRUNCATING PART
mutDATA <- read.csv (file = 'SOURCE_alltcga_jan2016.csv')

TRUNCATIONmutations = subset(mutDATA, mutation_type == 'Nonsense_Mutation') 

# first merge genes motif data with truncation mutations
merged <- merge(TRUNCATIONmutations, genbankkinases, by = 'Gene', nomatch = 0)


# will match ref amino acid as well as codon <endCatalytic) - AS A CHECKING MEASURE LOOK AT HOW THIS CHANGES IF JUST MATCHING BY CODON
# extract Codon / Ref  from amino acid change
merged$Codon <- as.numeric(str_extract(merged$amino_acid_change, "[0-9]+"))
merged$RefAmino <- str_extract(merged$amino_acid_change, "[A-Z]")


###############################################################################################################


# Locates the APE motif
numAPE <- str_count(APEfrag,  'APE')
genbankkinases$numAPE <- as.numeric(numAPE)
APELoc <- str_locate(APEfrag, 'APE')
APELoc <- as.data.frame(APELoc)
genbankkinases$APELoc <- as.numeric(APELoc$start)
genbankkinases$locAPEprotein <- genbankkinases$APEfragstart + genbankkinases$APELoc - 1

# Locates the PE motif
numPE <- str_count(APEfrag,  'PE')
genbankkinases$numPE <- as.numeric(numPE)
PELoc <- str_locate(APEfrag, 'PE')
PELoc <- as.data.frame(PELoc)
genbankkinases$PELoc <- as.numeric(PELoc$start)
genbankkinases$locPEprotein <- genbankkinases$APEfragstart + genbankkinases$PELoc - 2 # minus 2 as PE is 1 further back 



########################
# Determines correct APE and location information
########################

# Choose correct APE
genbankkinases$loc_actual_APE_A <- ifelse (genbankkinases$numAPE <1, genbankkinases$locGTx5NEprotein, genbankkinases$locAPEprotein)
genbankkinases$loc_actual_APE_B <- genbankkinases$loc_actual_APE_A
genbankkinases$loc_actual_APE_C   <- ifelse (genbankkinases$loc_actual_APE_A %in% NA, genbankkinases$locGTx6Eprotein, genbankkinases$loc_actual_APE_B)
genbankkinases$loc_actual_APE_D <- genbankkinases$loc_actual_APE_C 
genbankkinases$loc_actual_APE_E   <- ifelse (genbankkinases$loc_actual_APE_C %in% NA, genbankkinases$locPEprotein, genbankkinases$loc_actual_APE_D)
genbankkinases$loc_actual_APE_F <- genbankkinases$loc_actual_APE_E
genbankkinases$loc_actual_APE_G   <- ifelse (genbankkinases$loc_actual_APE_E %in% NA, genbankkinases$locGTx6Dprotein, genbankkinases$loc_actual_APE_F)
genbankkinases$loc_actual_APE_H <- genbankkinases$loc_actual_APE_G
genbankkinases$loc_actual_APE_I   <- ifelse (genbankkinases$loc_actual_APE_G %in% NA, genbankkinases$locGTxxYprotein, genbankkinases$loc_actual_APE_H)
genbankkinases$loc_actual_APE_J <- genbankkinases$loc_actual_APE_I
genbankkinases$loc_actual_APE_K   <- ifelse (genbankkinases$loc_actual_APE_I %in% NA, genbankkinases$locAPDprotein, genbankkinases$loc_actual_APE_J)
genbankkinases$loc_actual_APE_L  <- genbankkinases$loc_actual_APE_K
genbankkinases$loc_actual_APE_M   <- ifelse (genbankkinases$loc_actual_APE_K %in% NA, genbankkinases$locPPDprotein, genbankkinases$loc_actual_APE_L)
genbankkinases$loc_actual_APE_N  <- genbankkinases$loc_actual_APE_M
genbankkinases$loc_actual_APE_O   <- ifelse (genbankkinases$loc_actual_APE_M %in% NA, genbankkinases$WYxxPRLocprotein, genbankkinases$loc_actual_APE_N)
genbankkinases$loc_actual_APE_P <- genbankkinases$loc_actual_APE_O
genbankkinases$loc_actual_APE_Q   <- ifelse (genbankkinases$loc_actual_APE_O %in% NA, genbankkinases$locAxEprotein, genbankkinases$loc_actual_APE_P)
genbankkinases$loc_actual_APE_R <- genbankkinases$loc_actual_APE_Q
genbankkinases$loc_actual_APE_S   <- ifelse (genbankkinases$loc_actual_APE_Q %in% NA, genbankkinases$locPIRprotein, genbankkinases$loc_actual_APE_R)
genbankkinases$loc_actual_APE_T <- genbankkinases$loc_actual_APE_S
genbankkinases$loc_actual_APE_U   <- ifelse (genbankkinases$loc_actual_APE_S %in% NA, genbankkinases$locGx7Eprotein, genbankkinases$loc_actual_APE_T)
genbankkinases$loc_actual_APE_V <- genbankkinases$loc_actual_APE_U
genbankkinases$loc_actual_APE   <- ifelse (genbankkinases$loc_actual_APE_U %in% NA, genbankkinases$YxAPLocprotein, genbankkinases$loc_actual_APE_V)

genbankkinases$APE_A1loc   <- genbankkinases$loc_actual_APE
genbankkinases$APE_P1loc   <- genbankkinases$loc_actual_APE + 1
genbankkinases$APE_E1loc   <- genbankkinases$loc_actual_APE + 2
merged$endCAT <- as.numeric(merged$APE_E1loc)

###############################################################################################################


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
write.csv (freqSTRINGENTSHORTEST, file ='TEMP_freqTRANSCIPTSSHORTESTtcgaALL.csv')

#to length correct
freqSHORTEST <- read.csv ('TEMP_freqTRANSCIPTSSHORTESTtcgaALL.csv')
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

write.csv (freqADDEDSHORTEST, file ='TEMP_lengthcorrectedfreq_tcgaALLSHORTEST.csv')

#same for longest transcript
truncSTRINGENTLONGEST <- truncSTRINGENT[order(-truncSTRINGENT$endCAT),]
truncSTRINGENTLONGEST$duplicateREF <- paste(truncSTRINGENTLONGEST$Gene, truncSTRINGENTLONGEST$Codon, truncSTRINGENTLONGEST$case_id, sep ='_')
truncSTRINGENTLONGEST$mergedDUP= !duplicated(truncSTRINGENTLONGEST$duplicateREF)
genefreqLONGEST <- subset(truncSTRINGENTLONGEST, mergedDUP == 'TRUE')
genefreqLONGEST$transcriptFREQ <- paste(genefreqLONGEST$f, genefreqLONGEST$endCAT, sep = '_')
freqSTRINGENTLONGEST <-table(genefreqLONGEST$transcriptFREQ)
write.csv (freqSTRINGENTLONGEST, file ='TEMP_freqTRANSCIPTSLONGESTtcgaALL.csv')
freqLONGEST <- read.csv ('TEMP_freqTRANSCIPTSLONGESTtcgaALL.csv')
colnames(freqLONGEST)[2] <- "transcriptFREQ"
freqLONGEST$Gene <- gsub    ( "[-].*$"     , "", freqLONGEST$transcriptFREQ)
freqLONGEST$endCAT <-  gsub ( ".*[_]"    , "", freqLONGEST$transcriptFREQ)
freqLONGEST$endCAT <- as.numeric(freqLONGEST$endCAT)
freqLONGEST$oneoverendCAT <- 1 / freqLONGEST$endCAT
freqLONGEST$score <- freqLONGEST$Freq * freqLONGEST$oneoverendCAT
freqADDEDLONGEST <- ddply(freqLONGEST,"Gene",numcolwise(sum))
write.csv (freqADDEDLONGEST, file ='TEMP_lengthcorrectedfreq_tcgaALLLONGEST.csv')












# merge two SHORTEST truncating scores (using dplyr)
ccletruncSHORTEST <- read.csv (file = 'TEMP_lengthcorrectedfreq_ccleALLSHORTEST.csv')
tcgatruncSHORTEST <- read.csv (file = 'TEMP_lengthcorrectedfreq_tcgaALLSHORTEST.csv')
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
top30check <- rbind (validateTCGA, validateCCLE) 


#check the top30 to make sure the all transcript analysis doesnt give any false positives
top30 <- read.csv ('SOURCE_Top30.csv')
top30check <- merge (top30, top30check, by = 'Gene', nomatch = 0)
write.csv (top30check, file ='checking/top30check.csv')

#these were checked and are all valid - not extract catalytic frag of top 30 - note there are mutliple diff transcripts - use longest one that matches
# use genefreqLONGEST which is all tcga matches ranked by descending size of transcript
#first get rid of duplicates for gene to leave longest transcript
FRAGEXTRACT <- genefreqLONGEST
FRAGEXTRACT$DUP= !duplicated(FRAGEXTRACT$Gene)
FRAGEXTRACT <- subset(FRAGEXTRACT, DUP == 'TRUE')
FRAGEXTRACTtop30 <- merge (top30, FRAGEXTRACT, by = 'Gene', nomatch = 0)
#AMHR2 was changed manually as it is incorrect
write.csv (FRAGEXTRACTtop30, file ='TEMP_check.csv')

