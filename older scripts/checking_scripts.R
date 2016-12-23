# checking scripts

GENBANKKINASEScheck   <- read.csv ( file = 'SOURCE_genbankkinases.csv')
allTCGAcheck          <- read.csv ( file = 'newTCGA.csv')
CCLEcheck             <- read.csv ( file = 'SOURCE_ccletrimmed.csv')
COMBOcheck            <- read.csv ( file = 'SOURCE_combinedTCGAccle.csv')

# need to check the mutation input file - first check that that 433 kinases are listed in this data
  # start with TCGA (as CCLE may be annotated differently)


allTCGAcheck$DUP  <- !duplicated(allTCGAcheck$Gene)
uniqueTCGA   <- subset     (allTCGAcheck, DUP == 'TRUE')
uniqueTCGA   <- as.vector(uniqueTCGA$Gene)
write.csv (uniqueTCGA, file = 'outputs/checks/uniqueTCGA.csv')

genbanknames      <- GENBANKKINASEScheck
genbanknamesdup      <- GENBANKKINASEScheck
genbanknamesdup$DUP  <- !duplicated(genbanknames$Gene)
uniqueGENBANK   <- subset     (genbanknamesdup, DUP == 'TRUE')
uniqueGENBANK   <- as.vector(uniqueGENBANK$Gene)
write.csv (uniqueGENBANK, file = 'outputs/checks/uniqueGENBANK.csv')

TCGAvsGENBANK <- setdiff (uniqueGENBANK, uniqueTCGA)
write.csv (TCGAvsGENBANK, file = 'outputs/checks/TCGAvsGENBANK.csv')
### With the updated Cbio data and culls to ATM,ATR,SMG1 and MTOR these are missing, SBK3 is still missing but it looks like ryk now has some mutations
        # this shows that the only 2 of our 433 kinases not seen in the TCGA dataset is SBK3 and RYK
            # SBK3 does not appear to have any mutations in the TCGA dataset (according to cBio webportal)
            # ryk does have mutations - SO THIS DATA IS MISSING FROM THE TCGA WE LAST DOWNLOADED

  # now do the same with CCLE

CCLEcheck$DUP  <- !duplicated(CCLEcheck$Gene)
uniqueCCLE   <- subset     (CCLEcheck, DUP == 'TRUE')
uniqueCCLE   <- as.vector(uniqueCCLE$Gene)
write.csv (uniqueCCLE, file = 'outputs/checks/uniqueCCLE.csv')

CCLEvsGENBANK <- setdiff (uniqueGENBANK, uniqueCCLE)
write.csv (CCLEvsGENBANK, file = 'outputs/checks/CCLEvsGENBANK.csv')
    # map3k15 has no mutations
    # ysk should be recoded to map3k19
    # nim1k has no recorded mutations
    # mok has no mutations
    # neither has sbk3 nor sgk110 (the equivalent) are listed in ccle
    # mst4 should be changed to stk26
    # THESE HAVE NOW BEEN CHANGED in SOURCE_combinedTCGAccle.csv SO SHOULD NOT APPEAR AGAIN WHEN COMBO COMPARED BUT STILL PRESENT IN SOURCE_ccletrimmed so will appear here

COMBOcheck$DUP  <- !duplicated(COMBOcheck$Gene)
uniqueCOMBO   <- subset     (COMBOcheck, DUP == 'TRUE')
uniqueCOMBO   <- as.vector(uniqueCOMBO$Gene)
write.csv (uniqueCOMBO, file = 'outputs/checks/uniqueCOMBO.csv')

COMBOvsGENBANK <- setdiff (uniqueGENBANK, uniqueCOMBO)
write.csv (COMBOvsGENBANK, file = 'outputs/checks/COMBOvsGENBANK.csv')

# ryk disappears from combo as there must by ccle ryk mutations, likewise map3k15, nim1k, mok dissapears from combo as there must be tcga mutations
# sbk3 perisistently absence in both - THIS IS OK AS THERE IS NO MUTATIONAL DATA FOR SBK3 IN EITHER
## NEED TO RESOLVE RYK ISSUE IN TCGA AS MUTATIONAL DATA MISSING


# now need to check the correct transcript for each gene - so we can cut down on multiple mutations being reported from diff length transcripts
# because tcga and ccle transcripts may be different I think we will need to separate these (and not use a combofile for matching)
# start with TCGA

TCGAtranscipts <- allTCGAcheck

TCGAtranscipts$Codon <- str_extract(TCGAtranscipts$amino_acid_change, "[0-9]+")
TCGAtranscipts$Ref <- str_extract(TCGAtranscipts$amino_acid_change, "[A-Z]")

TCGAcrossref = merge(TCGAtranscipts, GENBANKKINASEScheck, by = 'Gene', nomatch = 0)

# cross ref codon with each protein seq of each respective transcript
codonCROSS <- as.numeric(TCGAcrossref$Codon)
TCGAcrossref$check <- str_sub(TCGAcrossref$Protein_Seq, codonCROSS, codonCROSS)
# now flag each one that matches
TCGAcrossref$flag <- ifelse (TCGAcrossref$check == TCGAcrossref$Ref, 'MATCH', 'NO' )
write.csv (TCGAcrossref, file = 'outputs/checks/TCGAcrossref.csv')

# now need to identify and remove (globally) any transcript with a no
TCGAcrossref$pastedMATCH <- paste(TCGAcrossref$f, TCGAcrossref$flag, sep = "_" )
TCGAcrossref$DUP= !duplicated(TCGAcrossref$pastedMATCH)
TCGAcrossrefNODUP <- subset(TCGAcrossref, DUP == 'TRUE')
# this gives answer for each gene id (f) if there is a MATCH and if there is a NO
# now want to select ids that are dont have a MATCH plus a NO (or vice versa) - to do this i think best is to seperate into 2 dataframes
TCGAcrossrefMATCH <- subset(TCGAcrossrefNODUP, flag == 'MATCH')
uniquematches <- TCGAcrossrefMATCH$Gene

uniquematches <- unique(uniquematches)
write.csv (uniquematches, file = 'outputs/checks/uniquematches.csv')
#### this gives 427 which appears to be 1 short - the missing ne is map3k14 which doesnt appear to have any missense mutations in the new dataset
# this gives 431 - which i think is correct given no muts in sbk3 and ryk issue

TCGAcrossrefNO <- subset(TCGAcrossrefNODUP, flag == 'NO')
write.csv (TCGAcrossrefNO, file = 'outputs/checks/TCGAcrossrefNO.csv')

matchedTCGAtranscripts <- as.vector(TCGAcrossrefMATCH$f)
unmatchedTCGAtranscripts <- as.vector(TCGAcrossrefNO$f)

TCGAnoNOs <- setdiff (matchedTCGAtranscripts, unmatchedTCGAtranscripts)
TCGAnoNOs <- as.data.frame (TCGAnoNOs)
colnames(TCGAnoNOs)[1] <- "f"

# so NoNos are a list of those transcripts with matches and no nos - there are multiple for some genes so need to choose the longest

TCGAnoNOs <- merge (TCGAnoNOs, GENBANKKINASEScheck , by = 'f', nomatch = 0)
colnames(TCGAnoNOs)[4] <- "Length"

# need to arrange by length
TCGAnoNOs <- TCGAnoNOs[order(-TCGAnoNOs$Length),]
# now look for duplicates - this should select the longest
TCGAnoNOs$DUP= !duplicated(TCGAnoNOs$Gene)
TCGAnoNOs <- subset(TCGAnoNOs, DUP == 'TRUE')

TCGAnoNOsnoDUPS <- as.data.frame(TCGAnoNOs$f)
colnames(TCGAnoNOsnoDUPS)[1] <- "f"
  
write.csv (TCGAnoNOsnoDUPS, file = 'outputs/checks/TCGAnoNOsnoDUPS.csv')

# now need to select these out of Source_GENBANK list which we use to do screen
GenbankTCGAcorrect <- merge (TCGAnoNOsnoDUPS, GENBANKKINASEScheck , by = 'f', nomatch = 0)
write.csv (GenbankTCGAcorrect, file = 'outputs/checks/GenbankTCGAcorrect.csv')

# this gives 118 which we know are correct - therefore need to fill in the others - the best way to do this will be to look at the matches and choose the transcript for each gene that is the most frequent and if there is a tie then choose the longest

TCGAmatches <- TCGAcrossrefMATCH$Gene
TCGAnoNOsGENES <- GenbankTCGAcorrect$Gene 

TCGAmatchesBUTnosTOO <- setdiff (TCGAmatches, TCGAnoNOsGENES)
TCGAmatchesBUTnosTOO <- as.data.frame (TCGAmatchesBUTnosTOO)
colnames(TCGAmatchesBUTnosTOO)[1] <- "Gene"

write.csv (TCGAmatchesBUTnosTOO, file = 'outputs/checks/TCGAmatchesBUTnosTOO.csv')
# this gives 90 (+341 = 431) which is the correct number minus sbk3 and ryk

TCGAcrossrefMATCHwithDUPS <- subset(TCGAcrossref, flag == 'MATCH')
TCGAmatchesBUTnoNONOs <- merge (TCGAcrossrefMATCHwithDUPS, TCGAmatchesBUTnosTOO , by = 'Gene', nomatch = 0)
write.csv (TCGAmatchesBUTnoNONOs, file = 'outputs/checks/TCGAmatchesBUTnoNONOs.csv')

uniqkinase <- TCGAmatchesBUTnoNONOs$f
unique(uniqkinase)


#need to remove any transcripts that do not have VAIK/HRD/DFG
# extract motifs
HRDmotifs  <- as.vector(TCGAmatchesBUTnoNONOs$HRDmotif)
VAIKmotifs <- as.vector(TCGAmatchesBUTnoNONOs$VAIKmotif)
DFGmotifs  <- as.vector(TCGAmatchesBUTnoNONOs$DFGmotif)
# count number of occurences of motif in each seq (using StringR package)
numberofDFG  <- str_count(TCGAmatchesBUTnoNONOs$Protein_Seq,  DFGmotifs)
numberofVAIK <- str_count(TCGAmatchesBUTnoNONOs$Protein_Seq, VAIKmotifs)
numberofHRD  <- str_count(TCGAmatchesBUTnoNONOs$Protein_Seq,  HRDmotifs)

TCGAmatchesBUTnoNONOs$numberofVAIK <- numberofVAIK
TCGAmatchesBUTnoNONOs$numberofHRD  <- numberofHRD
TCGAmatchesBUTnoNONOs$numberofDFG  <- numberofDFG

# remove any entry with number of VAIK or HRD or DFG = 0
TCGAmatchesBUTnoNONOsb4 <- TCGAmatchesBUTnoNONOs
write.csv (TCGAmatchesBUTnoNONOsb4, file = 'outputs/checks/TCGAmatchesBUTnoNONOsb4.csv')
uniqkinase <- TCGAmatchesBUTnoNONOsb4$Gene
unique(uniqkinase)


TCGAmatchesBUTnoNONOs <- subset(TCGAmatchesBUTnoNONOs, (numberofVAIK != 0))
TCGAmatchesBUTnoNONOs <- subset(TCGAmatchesBUTnoNONOs, (numberofHRD  != 0))
TCGAmatchesBUTnoNONOs <- subset(TCGAmatchesBUTnoNONOs, (numberofDFG  != 0))



freqOFmatches <-table(TCGAmatchesBUTnoNONOs$f)
write.csv (freqOFmatches, file = 'outputs/checks/freqOFmatches.csv')
# need to edit manually and then reload

TCGAfreqOfMATCHES          <- read.csv ( file = 'SOURCE_TCGAfreqOFmatchesEDIT2.csv')

# now remove duplicates and this should retain highest frequency
# first merge with something with the gene names

TCGAfreqOfMATCHESmerged <- merge (TCGAfreqOfMATCHES , TCGAcrossrefMATCHwithDUPS , by = 'f', nomatch = 0)

TCGAfreqOfMATCHESmerged <- subset(TCGAfreqOfMATCHESmerged, Freq != '0')
#this removes noNos from the list used to merge with

TCGAfreqOfMATCHESmerged$DUP= !duplicated(TCGAfreqOfMATCHESmerged$Gene)
TCGAfreqOfMATCHESmerged <- subset(TCGAfreqOfMATCHESmerged, DUP == 'TRUE')
write.csv (TCGAfreqOfMATCHESmerged, file = 'outputs/checks/TCGAfreqOfMATCHESmerged.csv')
# this gives 309 which is correct i think

# now just need to merge with the list of 118 to give master list of correct transcripts for TCGA

# then add all (all beginning with Methionine) transcripts of the missing 2 (SBK3 and RYK) incase further uses of TCGA dataset have mutations in these genes

# combine saved as SOURCE_Correct_TCGA_Transcript.csv - ACTUALLY NOW SOURCE_Correct_TCGA_Transcript2nd.csv





# now repeat for CCLE
GENBANKKINASEScheck   <- read.csv ( file = 'SOURCE_genbankkinases.csv')
CCLEtranscipts  <- read.csv ( file = 'SOURCE_ccletrimmed_namechanged.csv')

CCLEtranscipts$Codon <- str_extract(CCLEtranscipts$Protein_Change, "[0-9]+")
CCLEtranscipts$Ref <- str_extract(CCLEtranscipts$Protein_Change, "[A-Z]")

CCLEcrossref = merge(CCLEtranscipts, GENBANKKINASEScheck, by = 'Gene', nomatch = 0)

# cross ref codon with each protein seq of each respective transcript
codonCROSS <- as.numeric(CCLEcrossref$Codon)
CCLEcrossref$check <- str_sub(CCLEcrossref$Protein_Seq, codonCROSS, codonCROSS)
# now flag each one that matches
CCLEcrossref$flag <- ifelse (CCLEcrossref$check == CCLEcrossref$Ref, 'MATCH', 'NO' )
write.csv (CCLEcrossref, file = 'outputs/checks/CCLEcrossref.csv')

# now need to identify and remove (globally) any transcript with a no
CCLEcrossref$pastedMATCH <- paste(CCLEcrossref$f, CCLEcrossref$flag, sep = "_" )
CCLEcrossref$DUP= !duplicated(CCLEcrossref$pastedMATCH)
CCLEcrossrefNODUP <- subset(CCLEcrossref, DUP == 'TRUE')
# this gives answer for each gene id (f) if there is a MATCH and if there is a NO
# now want to select ids that are dont have a MATCH plus a NO (or vice versa) - to do this i think best is to seperate into 2 dataframes
CCLEcrossrefMATCH <- subset(CCLEcrossrefNODUP, flag == 'MATCH')
uniquematches <- CCLEcrossrefMATCH$Gene

uniquematches <- unique(uniquematches)
write.csv (uniquematches, file = 'outputs/checks/uniquematches.csv')
# this gives 429 - which i think is correct 

CCLEcrossrefNO <- subset(CCLEcrossrefNODUP, flag == 'NO')

matchedCCLEtranscripts <- as.vector(CCLEcrossrefMATCH$f)
unmatchedCCLEtranscripts <- as.vector(CCLEcrossrefNO$f)

CCLEnoNOs <- setdiff (matchedCCLEtranscripts, unmatchedCCLEtranscripts)
CCLEnoNOs <- as.data.frame (CCLEnoNOs)
colnames(CCLEnoNOs)[1] <- "f"

# so NoNos are a list of those transcripts with matches and no nos - there are multiple for some genes so need to choose the longest

CCLEnoNOs <- merge (CCLEnoNOs, GENBANKKINASEScheck , by = 'f', nomatch = 0)
colnames(CCLEnoNOs)[4] <- "Length"


# need to arrange by length
CCLEnoNOs <- CCLEnoNOs[order(-CCLEnoNOs$Length),]
# now look for duplicates - this should select the longest
CCLEnoNOs$DUP= !duplicated(CCLEnoNOs$Gene)
CCLEnoNOs <- subset(CCLEnoNOs, DUP == 'TRUE')

CCLEnoNOsnoDUPS <- as.data.frame(CCLEnoNOs$f)
colnames(CCLEnoNOsnoDUPS)[1] <- "f"

write.csv (CCLEnoNOsnoDUPS, file = 'outputs/checks/CCLEnoNOsnoDUPS.csv')

# now need to select these out of Source_GENBANK list which we use to do screen

GenbankCCLEcorrect <- merge (CCLEnoNOsnoDUPS, GENBANKKINASEScheck , by = 'f', nomatch = 0)
write.csv (GenbankCCLEcorrect, file = 'outputs/checks/GenbankCCLEcorrect.csv')

# this gives 395 which we know are correct - therefore need to fill in the others - the best way to do this will be to look at the matches and choose the transcript for each gene that is the most frequent and if there is a tie then choose the longest

CCLEmatches <- CCLEcrossrefMATCH$Gene
CCLEnoNOsGENES <- GenbankCCLEcorrect$Gene 

CCLEmatchesBUTnosTOO <- setdiff (CCLEmatches, CCLEnoNOsGENES)
CCLEmatchesBUTnosTOO <- as.data.frame (CCLEmatchesBUTnosTOO)
colnames(CCLEmatchesBUTnosTOO)[1] <- "Gene"

write.csv (CCLEmatchesBUTnosTOO, file = 'outputs/checks/CCLEmatchesBUTnosTOO.csv')
# this gives 34 (+395 = 429) which is the correct number minus the missing ones

CCLEcrossrefMATCHwithDUPS <- subset(CCLEcrossref, flag == 'MATCH')
CCLEmatchesBUTnoNONOs <- merge (CCLEcrossrefMATCHwithDUPS, CCLEmatchesBUTnosTOO , by = 'Gene', nomatch = 0)
write.csv (CCLEmatchesBUTnoNONOs, file = 'outputs/checks/CCLEmatchesBUTnoNONOs.csv')
freqOFmatches <-table(CCLEmatchesBUTnoNONOs$f)
write.csv (freqOFmatches, file = 'outputs/checks/freqOFmatches.csv')
# need to edit manually and then reload

CCLEfreqOfMATCHES          <- read.csv ( file = 'SOURCE_CCLEfreqOFmatchesEDIT.csv')

# now remove duplicates and this should retain highest frequency
# first merge with something with the gene names

CCLEfreqOfMATCHESmerged <- merge (CCLEfreqOfMATCHES , CCLEcrossrefMATCHwithDUPS , by = 'f', nomatch = 0)

CCLEfreqOfMATCHESmerged <- subset(CCLEfreqOfMATCHESmerged, Freq != '0')
#this removes noNos from the list used to merge with

CCLEfreqOfMATCHESmerged$DUP= !duplicated(CCLEfreqOfMATCHESmerged$Gene)
CCLEfreqOfMATCHESmerged <- subset(CCLEfreqOfMATCHESmerged, DUP == 'TRUE')
write.csv (CCLEfreqOfMATCHESmerged, file = 'outputs/checks/CCLEfreqOfMATCHESmerged.csv')
# this gives 34 which is the correct number

# now just need to merge with the list of 395 to give master list of correct transcripts for CCLE

##CCLE file saved as SOURCE_Correct_CCLE_Transcripts with column header f - this contains 429 - I havent added the missing ones as the CCLE dataset will not change - the TCGA one may change if we used data from a later date so will add the missing 2 (SBK3 and RYK)


