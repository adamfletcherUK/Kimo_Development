# TCGA script

#load StringR package
library("stringr", lib.loc="~/Library/R/3.1/library")
library("plyr", lib.loc="~/Library/R/3.1/library")

# match definitive gene names with GENBANK
genbank <- read.csv ( file = 'SOURCE_DEFINITIVE_GENBANK_TRIM.csv')

# merge with kinase gene names with motif data
names   <- read.csv (file = 'SOURCE_DEFINITIVEGenename_TRIM.csv') 
genbankkinases = merge(genbank, names, by = 'Gene', nomatch = 0)
write.csv (genbankkinases, 'check.csv')
#check how many unique kinases
uniqkinase <- genbankkinases$Gene
unique(uniqkinase)
# this displays 420 unique kinase names in TCGA

# extract motifs as per Manning et al
HRDmotifs  <- as.vector(genbankkinases$HRDmotif)
VAIKmotifs <- as.vector(genbankkinases$VAIKmotif)
DFGmotifs  <- as.vector(genbankkinases$DFGmotif)

# count number of occurences of motif in each seq (using StringR package)
numberofDFG  <- str_count(genbankkinases$Protein_Seq,  DFGmotifs)
numberofVAIK <- str_count(genbankkinases$Protein_Seq, VAIKmotifs)
numberofHRD  <- str_count(genbankkinases$Protein_Seq,  HRDmotifs)
genbankkinases$numberofVAIK <- numberofVAIK
genbankkinases$numberofHRD  <- numberofHRD
genbankkinases$numberofDFG  <- numberofDFG


# remove any entry with number of VAIK or HRD or DFG = 0
genbankkinases <- subset(genbankkinases, (numberofVAIK != 0))
genbankkinases <- subset(genbankkinases, (numberofHRD  != 0))
genbankkinases <- subset(genbankkinases, (numberofDFG  != 0))

# now need to visually inspect all kinases that have more than 1 of each motif, so that can call which one
# SOURCE_motif2use is a file with the correct motif to use when there are multiple
motifToUse     <- read.csv ('SOURCE_motif2use.csv')
genbankkinases <- merge    (genbankkinases, motifToUse, by = 'f')

#clean up genbankkinases
colnames(genbankkinases)[2] <- "Gene"
colnames(genbankkinases)[3] <- "Length"
colnames(genbankkinases)[4] <- "Protein_Seq"
colnames(genbankkinases)[6] <- "VAIKmotif"
colnames(genbankkinases)[7] <- "HRDmotif"
colnames(genbankkinases)[8] <- "DFGmotif"


#locate all motifs
HRDmotifs  <- as.vector(genbankkinases$HRDmotif)
VAIKmotifs <- as.vector(genbankkinases$VAIKmotif)
DFGmotifs  <- as.vector(genbankkinases$DFGmotif)

# don't need locate all for VAIK as there is only one gene with 2 vaik motifs (map4k1) and the valid one is the first
firstVAIK <- str_locate(genbankkinases$Protein_Seq, VAIKmotifs)
firstVAIK <- as.data.frame(firstVAIK)
genbankkinases$firstVAIK <- as.numeric(firstVAIK$start)
genbankkinases$actual_lysine <- as.numeric(genbankkinases$firstVAIK + 3)

#locate all HRD
firstHRD <- str_locate(genbankkinases$Protein_Seq, HRDmotifs)
firstHRD <- as.data.frame(firstHRD)
genbankkinases$firstHRD <- firstHRD$start
firstHRD <- as.numeric(firstHRD$start)
endofstring = as.numeric(firstHRD +1000000)

beyondfirstHRD <- substr(genbankkinases$Protein_Seq, firstHRD + 1, endofstring)

secondHRD <- str_locate(beyondfirstHRD, HRDmotifs)
secondHRD <- as.data.frame(secondHRD)
genbankkinases$secondHRD <- secondHRD$start
secondHRD = as.numeric(genbankkinases$secondHRD + firstHRD)
genbankkinases$secondHRD <- secondHRD

beyondsecondHRD <- substr(genbankkinases$Protein_Seq, secondHRD + 1, endofstring)

thirdHRD <- str_locate(beyondsecondHRD, HRDmotifs)
thirdHRD <- as.data.frame(thirdHRD)
genbankkinases$thirdHRD <- thirdHRD$start
thirdHRD = as.numeric(genbankkinases$thirdHRD + secondHRD)
genbankkinases$thirdHRD <- thirdHRD

#for DFG
firstDFG <- str_locate(genbankkinases$Protein_Seq, DFGmotifs)
firstDFG <- as.data.frame(firstDFG)
genbankkinases$firstDFG <- firstDFG$start
firstDFG <- as.numeric(firstDFG$start)
endofstring = as.numeric(firstDFG +1000000)

beyondfirstDFG <- substr(genbankkinases$Protein_Seq, firstDFG + 1, endofstring)

secondDFG <- str_locate(beyondfirstDFG, DFGmotifs)
secondDFG <- as.data.frame(secondDFG)
genbankkinases$secondDFG <- secondDFG$start
secondDFG = as.numeric(genbankkinases$secondDFG + firstDFG)
genbankkinases$secondDFG <- secondDFG

beyondsecondDFG <- substr(genbankkinases$Protein_Seq, secondDFG + 1, endofstring)

thirdDFG <- str_locate(beyondsecondDFG, DFGmotifs)
thirdDFG <- as.data.frame(thirdDFG)
genbankkinases$thirdDFG <- thirdDFG$start
thirdDFG = as.numeric(genbankkinases$thirdDFG + secondDFG)
genbankkinases$thirdDFG <- thirdDFG

# choose correct HRD
genbankkinases$actual_HRD_H <- ifelse (genbankkinases$HRDtouse == '1', genbankkinases$firstHRD, ifelse (genbankkinases$HRDtouse == '2', genbankkinases$secondHRD , genbankkinases$thirdHRD))
genbankkinases$actual_HRD_R <- genbankkinases$actual_HRD_H + 1
genbankkinases$actual_HRD_D <- genbankkinases$actual_HRD_H + 2


# choose correct DFG

genbankkinases$actual_DFG_D <- ifelse (genbankkinases$DFGtouse == '1', genbankkinases$firstDFG, ifelse (genbankkinases$DFGtouse == '2', genbankkinases$secondDFG , genbankkinases$thirdDFG))
genbankkinases$actual_DFG_F <- genbankkinases$actual_DFG_D + 1
genbankkinases$actual_DFG_G <- genbankkinases$actual_DFG_D + 2



# find APE motif
APEfragstart <- genbankkinases$actual_DFG_G + 10
genbankkinases$APEfragstart <- as.numeric(APEfragstart)
APEfragend <- genbankkinases$actual_DFG_G + 70
genbankkinases$APEfragend <- as.numeric(APEfragend)
APEfrag <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APEfragstart, genbankkinases$APEfragend )
genbankkinases$APEfrag <- APEfrag

# find APE
numAPE <- str_count(APEfrag,  'APE')
genbankkinases$numAPE <- as.numeric(numAPE)
APELoc <- str_locate(APEfrag, 'APE')
APELoc <- as.data.frame(APELoc)
genbankkinases$APELoc <- as.numeric(APELoc$start)
genbankkinases$locAPEprotein <- genbankkinases$APEfragstart + genbankkinases$APELoc - 1

# find PE
numPE <- str_count(APEfrag,  'PE')
genbankkinases$numPE <- as.numeric(numPE)
PELoc <- str_locate(APEfrag, 'PE')
PELoc <- as.data.frame(PELoc)
genbankkinases$PELoc <- as.numeric(PELoc$start)
genbankkinases$locPEprotein <- genbankkinases$APEfragstart + genbankkinases$PELoc - 2 # minus 2 as PE is 1 further back 

# find GTxxxxxxE
numGTx6E <- str_count(APEfrag,  'GT[A-Z]{6}E')
genbankkinases$numGTx6E  <- as.numeric(numGTx6E)
GTx6ELoc <- str_locate(APEfrag, 'GT[A-Z]{6}E')
GTx6ELoc <- as.data.frame(GTx6ELoc)
genbankkinases$GTx6ELoc <- as.numeric(GTx6ELoc$start)
genbankkinases$locGTx6Eprotein <- genbankkinases$APEfragstart + genbankkinases$GTx6ELoc + 5

# find GTxxxxxNE
numGTx5NE <- str_count(APEfrag,  'GT[A-Z]{5}NE')
genbankkinases$numGTx5NE  <- as.numeric(numGTx5NE)
GTx5NELoc <- str_locate(APEfrag, 'GT[A-Z]{5}NE')
GTx5NELoc <- as.data.frame(GTx5NELoc)
genbankkinases$GTx5NELoc <- as.numeric(GTx5NELoc$start)
genbankkinases$locGTx5NEprotein <- genbankkinases$APEfragstart + genbankkinases$GTx5NELoc + 5

# find GTxxxxxxD
numGTx6D <- str_count(APEfrag,  'GT[A-Z]{6}D')
genbankkinases$numGTx6D  <- as.numeric(numGTx6D)
GTx6DLoc <- str_locate(APEfrag, 'GT[A-Z]{6}D')
GTx6DLoc <- as.data.frame(GTx6DLoc)
genbankkinases$GTx6DLoc <- as.numeric(GTx6DLoc$start)
genbankkinases$locGTx6Dprotein <- genbankkinases$APEfragstart + genbankkinases$GTx6DLoc + 5

# find AxE
numAxE <- str_count(APEfrag,  'A[A-Z]E')
genbankkinases$numAxE <- as.numeric(numAxE)
AxELoc <- str_locate(APEfrag, 'A[A-Z]E')
AxELoc <- as.data.frame(AxELoc)
genbankkinases$AxELoc <- as.numeric(AxELoc$start)
genbankkinases$locAxEprotein <- genbankkinases$APEfragstart + genbankkinases$AxELoc - 1

# find APD
numAPD <- str_count(APEfrag,  'APD')
genbankkinases$numAPD  <- as.numeric(numAPD)
APDLoc <- str_locate(APEfrag, 'APD')
APDLoc <- as.data.frame(APDLoc)
genbankkinases$APDLoc <- as.numeric(APDLoc$start)
genbankkinases$locAPDprotein <- genbankkinases$APEfragstart + genbankkinases$APDLoc - 1

# find PPD
numPPD <- str_count(APEfrag,  'PPD')
genbankkinases$numPPD  <- as.numeric(numPPD)
PPDLoc <- str_locate(APEfrag, 'PPD')
PPDLoc <- as.data.frame(PPDLoc)
genbankkinases$PPDLoc <- as.numeric(PPDLoc$start)
genbankkinases$locPPDprotein <- genbankkinases$APEfragstart + genbankkinases$PPDLoc - 1

# find GTxxY
numGTxxY <- str_count(APEfrag,  'GT[A-Z]{2}Y')
genbankkinases$numGTxxY  <- as.numeric(numGTxxY)
GTxxYLoc <- str_locate(APEfrag, 'GT[A-Z]{2}Y')
GTxxYLoc <- as.data.frame(GTxxYLoc)
genbankkinases$GTxxYLoc <- as.numeric(GTxxYLoc$start)
genbankkinases$locGTxxYprotein <- genbankkinases$APEfragstart + genbankkinases$GTxxYLoc + 5

# find PIR
numPIR <- str_count(APEfrag,  'PIR')
genbankkinases$numPIR  <- as.numeric(numPIR)
PIRLoc <- str_locate(APEfrag, 'PIR')
PIRLoc <- as.data.frame(PIRLoc)
genbankkinases$PIRLoc <- as.numeric(PIRLoc$start)
genbankkinases$locPIRprotein <- genbankkinases$APEfragstart + genbankkinases$PIRLoc + 4

# find Gx7E
numGx7E <- str_count(APEfrag,  'G[A-Z]{7}E')
genbankkinases$numGx7E  <- as.numeric(numGx7E)
Gx7ELoc <- str_locate(APEfrag, 'G[A-Z]{7}E')
Gx7ELoc <- as.data.frame(Gx7ELoc)
genbankkinases$Gx7ELoc <- as.numeric(Gx7ELoc$start)
genbankkinases$locGx7Eprotein <- genbankkinases$APEfragstart + genbankkinases$Gx7ELoc + 5

# find YxAP (captures MAPKAPK5)
numYxAP <- str_count(APEfrag,  'Y[A-Z]{1}AP')
genbankkinases$numYxAP  <- as.numeric(numYxAP)
YxAPLoc <- str_locate(APEfrag, 'Y[A-Z]{1}AP')
YxAPLoc <- as.data.frame(YxAPLoc)
genbankkinases$YxAPLoc <- as.numeric(YxAPLoc$start)
genbankkinases$YxAPLocprotein <- genbankkinases$APEfragstart + genbankkinases$YxAPLoc + 1

# find WYxxPR (captures MAPK4 and MAPK6)
numWYxxPR <- str_count(APEfrag,  'WY[A-Z]{2}PR')
genbankkinases$numWYxxPR  <- as.numeric(numWYxxPR)
WYxxPRLoc <- str_locate(APEfrag, 'WY[A-Z]{2}PR')
WYxxPRLoc <- as.data.frame(WYxxPRLoc)
genbankkinases$WYxxPRLoc <- as.numeric(WYxxPRLoc$start)
genbankkinases$WYxxPRLocprotein <- genbankkinases$APEfragstart + genbankkinases$WYxxPRLoc + 2


# choose correct APE
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


genbankkinases$APE_m14loc   <- genbankkinases$loc_actual_APE - 14
genbankkinases$APE_m13loc   <- genbankkinases$loc_actual_APE - 13
genbankkinases$APE_m12loc   <- genbankkinases$loc_actual_APE - 12
genbankkinases$APE_m11loc   <- genbankkinases$loc_actual_APE - 11
genbankkinases$APE_m10loc   <- genbankkinases$loc_actual_APE - 10
genbankkinases$APE_m9loc   <- genbankkinases$loc_actual_APE - 9
genbankkinases$APE_m8loc   <- genbankkinases$loc_actual_APE - 8
genbankkinases$APE_m7loc   <- genbankkinases$loc_actual_APE - 7
genbankkinases$APE_m6loc   <- genbankkinases$loc_actual_APE - 6
genbankkinases$APE_m5loc   <- genbankkinases$loc_actual_APE - 5
genbankkinases$APE_m4loc   <- genbankkinases$loc_actual_APE - 4
genbankkinases$APE_m3loc   <- genbankkinases$loc_actual_APE - 3
genbankkinases$APE_m2loc   <- genbankkinases$loc_actual_APE - 2
genbankkinases$APE_m1loc   <- genbankkinases$loc_actual_APE - 1
genbankkinases$APE_A1loc   <- genbankkinases$loc_actual_APE
genbankkinases$APE_P1loc   <- genbankkinases$loc_actual_APE + 1
genbankkinases$APE_E1loc   <- genbankkinases$loc_actual_APE + 2




genbankkinases$APE_m14       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m14loc , genbankkinases$APE_m14loc)
genbankkinases$APE_m13       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m13loc , genbankkinases$APE_m13loc)
genbankkinases$APE_m12       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m12loc , genbankkinases$APE_m12loc)
genbankkinases$APE_m11       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m11loc , genbankkinases$APE_m11loc)
genbankkinases$APE_m10       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m10loc , genbankkinases$APE_m10loc)
genbankkinases$APE_m9       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m9loc , genbankkinases$APE_m9loc)
genbankkinases$APE_m8       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m8loc , genbankkinases$APE_m8loc)
genbankkinases$APE_m7       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m7loc , genbankkinases$APE_m7loc)
genbankkinases$APE_m6       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m6loc , genbankkinases$APE_m6loc)
genbankkinases$APE_m5       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m5loc , genbankkinases$APE_m5loc)
genbankkinases$APE_m4       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m4loc , genbankkinases$APE_m4loc)
genbankkinases$APE_m3       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m3loc , genbankkinases$APE_m3loc)
genbankkinases$APE_m2       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m2loc , genbankkinases$APE_m2loc)
genbankkinases$APE_m1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m1loc , genbankkinases$APE_m1loc)
genbankkinases$APE_A1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_A1loc , genbankkinases$APE_A1loc)
genbankkinases$APE_P1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_P1loc , genbankkinases$APE_P1loc)
genbankkinases$APE_E1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_E1loc , genbankkinases$APE_E1loc)

genbankkinases$APE_MOTIF <- paste(genbankkinases$APE_A1, genbankkinases$APE_P1, genbankkinases$APE_E1)
genbankkinases$PE_MOTIF <- paste(genbankkinases$APE_P1, genbankkinases$APE_E1)

# if no G in APEm6 then look for G in APEm8 and then APEm9 and then A's in these positions

genbankkinases$correctedG_loc <- ifelse (genbankkinases$APE_m6 != 'G', genbankkinases$APE_m9loc, genbankkinases$APE_m6loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG != 'G', genbankkinases$APE_m5loc, genbankkinases$correctedG_loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG != 'G', genbankkinases$APE_m8loc, genbankkinases$correctedG_loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG != 'G', genbankkinases$APE_m10loc, genbankkinases$correctedG_loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG != 'G', genbankkinases$APE_m11loc, genbankkinases$correctedG_loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG != 'G', genbankkinases$APE_m7loc, genbankkinases$correctedG_loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG == 'G', genbankkinases$correctedG_loc, genbankkinases$APE_m6loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

#end of kinase domain = APE + defined end
genbankkinases$endCAT <- as.numeric(genbankkinases$APE_E1loc)

# find GxGxxG
GxGfragstart <- genbankkinases$actual_lysine - 45
genbankkinases$GxGfragstart <- as.numeric(GxGfragstart)

# NEED TO PUT THIS LINE IN FOR ALL MOTIF SEARCHES - TELLING THE FRAG TO START AT 1 IF MINUS NUMBER
genbankkinases$GxGfragstart <- ifelse(genbankkinases$GxGfragstart < 1, 1,genbankkinases$GxGfragstart )

GxGfragend <- genbankkinases$actual_lysine - 4
genbankkinases$GxGfragend <- as.numeric(GxGfragend)
GxGfrag <- str_sub(genbankkinases$Protein_Seq, genbankkinases$GxGfragstart, genbankkinases$GxGfragend )
genbankkinases$GxGfrag <- GxGfrag

# find GxGxxG
numGxGxxG <- str_count(GxGfrag,  'G[A-Z]G[A-Z]{2}G')
genbankkinases$numGxGxxG <- as.numeric(numGxGxxG)
GxGxxGLoc <- str_locate(GxGfrag, 'G[A-Z]G[A-Z]{2}G')
GxGxxGLoc <- as.data.frame(GxGxxGLoc)
genbankkinases$GxGxxGLoc <- as.numeric(GxGxxGLoc$start)
genbankkinases$locGxGxxGprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGxxGLoc - 1

# find GxGxF
numGxGxF <- str_count(GxGfrag,  'G[A-Z]G[A-Z]F')
genbankkinases$numGxGxF <- as.numeric(numGxGxF)
GxGxFLoc <- str_locate(GxGfrag, 'G[A-Z]G[A-Z]F')
GxGxFLoc <- as.data.frame(GxGxFLoc)
genbankkinases$GxGxFLoc <- as.numeric(GxGxFLoc$start)
genbankkinases$locGxGxFprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGxFLoc - 1

# find GxGxY
numGxGxY <- str_count(GxGfrag,  'G[A-Z]G[A-Z]Y')
genbankkinases$numGxGxY <- as.numeric(numGxGxY)
GxGxYLoc <- str_locate(GxGfrag, 'G[A-Z]G[A-Z]Y')
GxGxYLoc <- as.data.frame(GxGxYLoc)
genbankkinases$GxGxYLoc <- as.numeric(GxGxYLoc$start)
genbankkinases$locGxGxYprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGxYLoc - 1

# find GxFG
numGxFG <- str_count(GxGfrag,  'G[A-Z]FG')
genbankkinases$numGxFG <- as.numeric(numGxFG)
GxFGLoc <- str_locate(GxGfrag, 'G[A-Z]FG')
GxFGLoc <- as.data.frame(GxFGLoc)
genbankkinases$GxFGLoc <- as.numeric(GxFGLoc$start)
genbankkinases$locGxFGprotein <- genbankkinases$GxGfragstart + genbankkinases$GxFGLoc - 3
#above is minus 3 rather than minus 1 to take account of the GxFG starting 2 on from GxGxxG

# find GxGxxA
numGxGxxA <- str_count(GxGfrag,  'G[A-Z]G[A-Z]{2}A')
genbankkinases$numGxGxxA <- as.numeric(numGxGxxA)
GxGxxALoc <- str_locate(GxGfrag, 'G[A-Z]G[A-Z]{2}A')
GxGxxALoc <- as.data.frame(GxGxxALoc)
genbankkinases$GxGxxALoc <- as.numeric(GxGxxALoc$start)
genbankkinases$locGxGxxAprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGxxALoc - 1

# find GxGxxS
numGxGxxS <- str_count(GxGfrag,  'G[A-Z]G[A-Z]{2}S')
genbankkinases$numGxGxxS <- as.numeric(numGxGxxS)
GxGxxSLoc <- str_locate(GxGfrag, 'G[A-Z]G[A-Z]{2}S')
GxGxxSLoc <- as.data.frame(GxGxxSLoc)
genbankkinases$GxGxxSLoc <- as.numeric(GxGxxSLoc$start)
genbankkinases$locGxGxxSprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGxxSLoc - 1

# find SxGxxG
numSxGxxG <- str_count(GxGfrag,  'S[A-Z]G[A-Z]{2}G')
genbankkinases$numSxGxxG <- as.numeric(numSxGxxG)
SxGxxGLoc <- str_locate(GxGfrag, 'S[A-Z]G[A-Z]{2}G')
SxGxxGLoc <- as.data.frame(SxGxxGLoc)
genbankkinases$SxGxxGLoc <- as.numeric(SxGxxGLoc$start)
genbankkinases$locSxGxxGprotein <- genbankkinases$GxGfragstart + genbankkinases$SxGxxGLoc - 1

# find GxxG
numGxxG <- str_count(GxGfrag,  'G[A-Z]{2}G')
genbankkinases$numGxxG <- as.numeric(numGxxG)
GxxGLoc <- str_locate(GxGfrag, 'G[A-Z]{2}G')
GxxGLoc <- as.data.frame(GxxGLoc)
genbankkinases$GxxGLoc <- as.numeric(GxxGLoc$start)
genbankkinases$locGxxGprotein <- genbankkinases$GxGfragstart + genbankkinases$GxxGLoc - 3

# find G4xG
numG4xG <- str_count(GxGfrag,  'G[A-Z]{4}G')
genbankkinases$numG4xG <- as.numeric(numG4xG)
G4xGLoc <- str_locate(GxGfrag, 'G[A-Z]{4}G')
G4xGLoc <- as.data.frame(G4xGLoc)
genbankkinases$G4xGLoc <- as.numeric(G4xGLoc$start)
genbankkinases$locG4xGprotein <- genbankkinases$GxGfragstart + genbankkinases$G4xGLoc - 1

# find GxG
numGxG <- str_count(GxGfrag,  'G[A-Z]G')
genbankkinases$numGxG <- as.numeric(numGxG)
GxGLoc <- str_locate(GxGfrag, 'G[A-Z]G')
GxGLoc <- as.data.frame(GxGLoc)
genbankkinases$GxGLoc <- as.numeric(GxGLoc$start)
genbankkinases$locGxGprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGLoc - 1

# find GxTxF (locates PIK3R4 - need to check it is the correct motif by looking at structure)
numGxTxF <- str_count(GxGfrag,  'G[A-Z]T[A-Z]F')
genbankkinases$numGxTxF <- as.numeric(numGxTxF)
GxTxFLoc <- str_locate(GxGfrag, 'G[A-Z]T[A-Z]F')
GxTxFLoc <- as.data.frame(GxTxFLoc)
genbankkinases$GxTxFLoc <- as.numeric(GxTxFLoc$start)
genbankkinases$locGxTxFprotein <- genbankkinases$GxGfragstart + genbankkinases$GxTxFLoc - 1

# find GG (only use as the last resort to capture AMHRH2 and ANKK1)
numGG <- str_count(GxGfrag, 'GG')
genbankkinases$numGG <- as.numeric(numGG)
GGLoc <- str_locate(GxGfrag, 'GG')
GGLoc <- as.data.frame(GGLoc)
genbankkinases$GGLoc <- as.numeric(GGLoc$start)
genbankkinases$locGGprotein <- genbankkinases$GxGfragstart + genbankkinases$GGLoc - 3


# choose correct GxG
genbankkinases$loc_actual_gxg_A <- ifelse (genbankkinases$numGxGxxG <1, genbankkinases$locGxGxFprotein, genbankkinases$locGxGxxGprotein)
genbankkinases$loc_actual_gxg_B <- genbankkinases$loc_actual_gxg_A
genbankkinases$loc_actual_gxg_C   <- ifelse (genbankkinases$loc_actual_gxg_A %in% NA, genbankkinases$locGxFGprotein, genbankkinases$loc_actual_gxg_B)
genbankkinases$loc_actual_gxg_D <- genbankkinases$loc_actual_gxg_C
genbankkinases$loc_actual_gxg_E   <- ifelse (genbankkinases$loc_actual_gxg_C %in% NA, genbankkinases$locGxGxYprotein, genbankkinases$loc_actual_gxg_D)
genbankkinases$loc_actual_gxg_F <- genbankkinases$loc_actual_gxg_E

genbankkinases$loc_actual_gxg_G   <- ifelse (genbankkinases$loc_actual_gxg_E %in% NA, genbankkinases$locGxGxxAprotein, genbankkinases$loc_actual_gxg_F)
genbankkinases$loc_actual_gxg_H <- genbankkinases$loc_actual_gxg_G
genbankkinases$loc_actual_gxg_I   <- ifelse (genbankkinases$loc_actual_gxg_G %in% NA, genbankkinases$locGxGxxSprotein, genbankkinases$loc_actual_gxg_H)
genbankkinases$loc_actual_gxg_J <- genbankkinases$loc_actual_gxg_I
genbankkinases$loc_actual_gxg_K   <- ifelse (genbankkinases$loc_actual_gxg_I %in% NA, genbankkinases$locSxGxxGprotein, genbankkinases$loc_actual_gxg_J)
genbankkinases$loc_actual_gxg_L <- genbankkinases$loc_actual_gxg_K
genbankkinases$loc_actual_gxg_M   <- ifelse (genbankkinases$loc_actual_gxg_K %in% NA, genbankkinases$locGxxGprotein, genbankkinases$loc_actual_gxg_L)
genbankkinases$loc_actual_gxg_N <- genbankkinases$loc_actual_gxg_M
genbankkinases$loc_actual_gxg_O   <- ifelse (genbankkinases$loc_actual_gxg_M %in% NA, genbankkinases$locG4xGprotein, genbankkinases$loc_actual_gxg_N)
genbankkinases$loc_actual_gxg_P <- genbankkinases$loc_actual_gxg_O
genbankkinases$loc_actual_gxg_Q   <- ifelse (genbankkinases$loc_actual_gxg_O %in% NA, genbankkinases$locGxGprotein, genbankkinases$loc_actual_gxg_P)
genbankkinases$loc_actual_gxg_R <- genbankkinases$loc_actual_gxg_Q
genbankkinases$loc_actual_gxg_S   <- ifelse (genbankkinases$loc_actual_gxg_Q %in% NA, genbankkinases$locGxTxFprotein, genbankkinases$loc_actual_gxg_R)
genbankkinases$loc_actual_gxg_T <- genbankkinases$loc_actual_gxg_S
genbankkinases$loc_actual_gxg   <- ifelse (genbankkinases$loc_actual_gxg_S %in% NA, genbankkinases$locGGprotein, genbankkinases$loc_actual_gxg_T)

#extract GxGxxG motif and location/aa for each
startGxG <- genbankkinases$loc_actual_gxg
endGxG <- genbankkinases$loc_actual_gxg + 5
genbankkinases$GxGxxGmotif <- str_sub(genbankkinases$Protein_Seq, startGxG, endGxG)

genbankkinases$GxGxxG_G1     <- str_sub(genbankkinases$GxGxxGmotif, 1, 1)
genbankkinases$GxGxxG_X1     <- str_sub(genbankkinases$GxGxxGmotif, 2, 2)
genbankkinases$GxGxxG_G2     <- str_sub(genbankkinases$GxGxxGmotif, 3, 3)
genbankkinases$GxGxxG_X2     <- str_sub(genbankkinases$GxGxxGmotif, 4, 4)
genbankkinases$GxGxxG_X3     <- str_sub(genbankkinases$GxGxxGmotif, 5, 5)
genbankkinases$GxGxxG_G3     <- str_sub(genbankkinases$GxGxxGmotif, 6, 6)


genbankkinases$GxGxxG_G1loc   <- genbankkinases$loc_actual_gxg
genbankkinases$GxGxxG_X1loc   <- genbankkinases$loc_actual_gxg + 1
genbankkinases$GxGxxG_G2loc   <- genbankkinases$loc_actual_gxg + 2
genbankkinases$GxGxxG_X2loc   <- genbankkinases$loc_actual_gxg + 3
genbankkinases$GxGxxG_X3loc   <- genbankkinases$loc_actual_gxg + 4
genbankkinases$GxGxxG_G3loc   <- genbankkinases$loc_actual_gxg + 5



genbankkinases$GXGXXG_MOTIF <- paste(genbankkinases$GxGxxG_G1, genbankkinases$GxGxxG_G2, genbankkinases$GxGxxG_G3)

genbankkinases$CATALYTIC_FRAG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$GxGxxG_G1loc, genbankkinases$APE_E1loc )

write.csv (genbankkinases, file ='outputs/TCGAgenbankkinases.csv')

#find salt bridge e
genbankkinases$saltbridgeEfragSTART <- genbankkinases$actual_lysine + 5
genbankkinases$saltbridgeEfragEND <- genbankkinases$actual_lysine + 30
genbankkinases$saltbridgeEfrag <- str_sub(genbankkinases$Protein_Seq, genbankkinases$saltbridgeEfragSTART, genbankkinases$saltbridgeEfragEND )
saltbridgeEfrag <- genbankkinases$saltbridgeEfrag

#find ExxxL
numExxxL <- str_count(saltbridgeEfrag,  'E[A-Z]{3}L')
genbankkinases$numExxxL <- as.numeric(numExxxL)
ExxxLLoc <- str_locate(saltbridgeEfrag,  'E[A-Z]{3}L')
ExxxLLoc <- as.data.frame(ExxxLLoc)
genbankkinases$ExxxLLoc <- as.numeric(ExxxLLoc$start)
genbankkinases$ExxxLLocprotein <- genbankkinases$saltbridgeEfragSTART + genbankkinases$ExxxLLoc - 1

#find ExxI
numExxI <- str_count(saltbridgeEfrag,  'E[A-Z]{2}I')
genbankkinases$numExxI <- as.numeric(numExxI)
ExxILoc <- str_locate(saltbridgeEfrag,  'E[A-Z]{2}I')
ExxILoc <- as.data.frame(ExxILoc)
genbankkinases$ExxILoc <- as.numeric(ExxILoc$start)
genbankkinases$ExxILocprotein <- genbankkinases$saltbridgeEfragSTART + genbankkinases$ExxILoc - 1

#find ExxxM
numExxxM <- str_count(saltbridgeEfrag,  'E[A-Z]{3}M')
genbankkinases$numExxxM <- as.numeric(numExxxM)
ExxxMLoc <- str_locate(saltbridgeEfrag,  'E[A-Z]{3}M')
ExxxMLoc <- as.data.frame(ExxxMLoc)
genbankkinases$ExxxMLoc <- as.numeric(ExxxMLoc$start)
genbankkinases$ExxxMLocprotein <- genbankkinases$saltbridgeEfragSTART + genbankkinases$ExxxMLoc - 1

#find ExxL
numExxL <- str_count(saltbridgeEfrag,  'E[A-Z]{2}L')
genbankkinases$numExxL <- as.numeric(numExxL)
ExxLLoc <- str_locate(saltbridgeEfrag,  'E[A-Z]{2}L')
ExxLLoc <- as.data.frame(ExxLLoc)
genbankkinases$ExxLLoc <- as.numeric(ExxLLoc$start)
genbankkinases$ExxLLocprotein <- genbankkinases$saltbridgeEfragSTART + genbankkinases$ExxLLoc - 1

#find QxxxE
numQxxxE <- str_count(saltbridgeEfrag,  'Q[A-Z]{3}E')
genbankkinases$numQxxxE <- as.numeric(numQxxxE)
QxxxELoc <- str_locate(saltbridgeEfrag,  'Q[A-Z]{3}E')
QxxxELoc <- as.data.frame(QxxxELoc)
genbankkinases$QxxxELoc <- as.numeric(QxxxELoc$start)
genbankkinases$QxxxELocprotein <- genbankkinases$saltbridgeEfragSTART + genbankkinases$QxxxELoc + 3

# choose correct saltbridge
genbankkinases$loc_actual_saltbridge_A <- ifelse (genbankkinases$numExxxL <1, genbankkinases$ExxILocprotein, genbankkinases$ExxxLLocprotein)
genbankkinases$loc_actual_saltbridge_B <- genbankkinases$loc_actual_saltbridge_A
genbankkinases$loc_actual_saltbridge_C   <- ifelse (genbankkinases$loc_actual_saltbridge_A %in% NA, genbankkinases$ExxxMLocprotein, genbankkinases$loc_actual_saltbridge_B)
genbankkinases$loc_actual_saltbridge_D <- genbankkinases$loc_actual_saltbridge_C
genbankkinases$loc_actual_saltbridge_E   <- ifelse (genbankkinases$loc_actual_saltbridge_C %in% NA, genbankkinases$ExxLLocprotein, genbankkinases$loc_actual_saltbridge_D)
genbankkinases$loc_actual_saltbridge_F <- genbankkinases$loc_actual_saltbridge_E
genbankkinases$loc_actual_saltbridge   <- ifelse (genbankkinases$loc_actual_saltbridge_E %in% NA, genbankkinases$QxxxELocprotein, genbankkinases$loc_actual_saltbridge_F)

genbankkinases$checkE <- str_sub(genbankkinases$Protein_Seq, genbankkinases$loc_actual_saltbridge, genbankkinases$loc_actual_saltbridge)

write.csv (genbankkinases, file ='TEMP_genbankkinases.csv')

# merge with combined TCGA data
comboDATA <- read.csv (file = 'SOURCE_allTCGA_jan2016.csv')  #was'SOURCE_allTCGA.csv'# and can be allcbioNov2015.csv
kinaseLoc <- genbankkinases


comboDATA = subset(comboDATA, mutation_type != '3UTR')
comboDATA = subset(comboDATA, mutation_type != '5UTR')
comboDATA = subset(comboDATA, mutation_type != 'Intron')
comboDATA = subset(comboDATA, mutation_type != 'Splice_Site_Ins')
comboDATA = subset(comboDATA, mutation_type != 'Splice_Site_Del')
comboDATA = subset(comboDATA, mutation_type != 'Splice_Site_SNP')
comboDATA = subset(comboDATA, mutation_type != 'Silent')
comboDATA = subset(comboDATA, mutation_type != 'Nonsense_Mutation')
comboDATA = subset(comboDATA, mutation_type != 'Frame_Shift_Del')
comboDATA = subset(comboDATA, mutation_type != 'Frame_Shift_Ins')
comboDATA = subset(comboDATA, mutation_type != 'In_Frame_Del')
comboDATA = subset(comboDATA, mutation_type != 'In_Frame_Ins')
comboDATA = subset(comboDATA, mutation_type != 'outofframe')
comboDATA = subset(comboDATA, mutation_type != 'Flank')
comboDATA = subset(comboDATA, mutation_type != 'Splice_Site')
comboDATA = subset(comboDATA, mutation_type != 'De_novo_Start_InFrame')
comboDATA = subset(comboDATA, mutation_type != 'Stop_Codon_DNP')
comboDATA = subset(comboDATA, mutation_type != 'Stop_Codon_Ins')


#extract data from comboDATA
CodonCombo <- str_extract(comboDATA$amino_acid_change, "[0-9]+")
comboDATA$Codon <- CodonCombo
RefAmino <- str_extract(comboDATA$amino_acid_change, "[A-Z]")
comboDATA$AA <- RefAmino
VarAmino <- str_sub(comboDATA$amino_acid_change, -1, -1)
comboDATA$VarAA <- VarAmino


#GXGXXG_G1
kinaseLoc$Codon <- kinaseLoc$GxGxxG_G1loc
kinaseLoc$AA <- kinaseLoc$GxGxxG_G1
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_G1 <- merged
write.csv (merged_G1, file='outputs/tcga/merged_GxGxxG_G1_all.csv')

#GXGXXG_G2
kinaseLoc$Codon <- kinaseLoc$GxGxxG_G2loc
kinaseLoc$AA <- kinaseLoc$GxGxxG_G2
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_G2 <- merged
write.csv (merged_G2, file='outputs/tcga/merged_GxGxxG_G2_all.csv')

# critical lysine
kinaseLoc$Codon <- kinaseLoc$actual_lysine
kinaseLoc$AA <- 'K'
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_K <- merged
write.csv (merged_K, file='outputs/tcga/critical_lysine_all.csv')



# DFG_D
kinaseLoc$Codon <- kinaseLoc$actual_DFG_D
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_DFG_D, kinaseLoc$actual_DFG_D)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_DFG_D <- merged
write.csv (merged_DFG_D, file='outputs/tcga/merged_DFG_G_all.csv')

# DFG_F
kinaseLoc$Codon <- kinaseLoc$actual_DFG_F
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_DFG_F, kinaseLoc$actual_DFG_F)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_DFG_F <- merged
write.csv (merged_DFG_D, file='outputs/tcga/merged_DFG_F_all.csv')

# DFG_G
kinaseLoc$Codon <- kinaseLoc$actual_DFG_G
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_DFG_G, kinaseLoc$actual_DFG_G)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_DFG_G <- merged
write.csv (merged_DFG_G, file='outputs/tcga/merged_DFG_G_all.csv')


# HRD_H
kinaseLoc$Codon <- kinaseLoc$actual_HRD_H
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_HRD_H, kinaseLoc$actual_HRD_H)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_HRD_H <- merged
write.csv (merged_HRD_H, file='outputs/tcga/merged_HRD_H_all.csv')

# HRD_R
kinaseLoc$Codon <- kinaseLoc$actual_HRD_R
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_HRD_R, kinaseLoc$actual_HRD_R)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_HRD_R <- merged
write.csv (merged_HRD_R, file='outputs/tcga/merged_HRD_R_all.csv')

# HRD_D
kinaseLoc$Codon <- kinaseLoc$actual_HRD_D
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_HRD_D, kinaseLoc$actual_HRD_D)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_HRD_D <- merged
write.csv (merged_HRD_D, file='outputs/tcga/merged_HRD_D_all.csv')

# HRD + 5
kinaseLoc$Codon <- kinaseLoc$actual_HRD_D + 5
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_HRD_D, kinaseLoc$actual_HRD_D)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_HRD_Dp5 <- merged
write.csv (merged_HRD_Dp5, file='outputs/tcga/merged_HRD_D_all.csv')


# APE_A
kinaseLoc$Codon <- kinaseLoc$APE_A1loc 
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$APE_A1loc, kinaseLoc$APE_A1loc)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_APE_A <- merged
write.csv (merged_APE_A, file='outputs/tcga/merged_APE_A_all.csv')

# APE_P
kinaseLoc$Codon <- kinaseLoc$APE_P1loc 
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$APE_P1loc, kinaseLoc$APE_P1loc)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_APE_P <- merged
write.csv (merged_APE_P, file='outputs/tcga/merged_APE_P_all.csv')

# APE_E
kinaseLoc$Codon <- kinaseLoc$APE_E1loc 
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$APE_E1loc, kinaseLoc$APE_E1loc)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_APE_E <- merged
write.csv (merged_APE_E, file='outputs/tcga/merged_APE_E_all.csv')


# CorrectedG
kinaseLoc$Codon <- kinaseLoc$correctedG_loc
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$correctedG_loc, kinaseLoc$correctedG_loc)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_correctedG <- merged
write.csv (merged_correctedG, file='outputs/tcga/merged_correctedG_all.csv')

# SALTBRIDGE_E
kinaseLoc$Codon <- kinaseLoc$loc_actual_saltbridge
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$loc_actual_saltbridge, kinaseLoc$loc_actual_saltbridge)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_saltbridge <- merged
write.csv (merged_saltbridge, file='outputs/tcga/merged_saltbridge_all.csv')

# critical lysine minus 2 = V[A]IK
kinaseLoc$Codon <- kinaseLoc$actual_lysine - 2
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$Codon, kinaseLoc$Codon)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_Kminus2 <- merged
write.csv (merged_Kminus2, file='outputs/tcga/critical_lysine_minus2.csv')

#bind tables together
total <- rbind (merged_G2, merged_Kminus2, merged_HRD_H, merged_HRD_R, merged_HRD_D, merged_HRD_Dp5, merged_DFG_D, merged_DFG_G, merged_APE_A, merged_APE_P, merged_APE_E, merged_saltbridge, merged_correctedG)
write.csv (total, file ='outputs/tcga/total_ALLTRANSCRIPTS.csv')


freqscreen <- table(total$Gene)
write.csv (freqscreen, file ='SOURCE_tcgafreqscreen_ALLTRANSCRIPTS.csv')         

#merge CCLE and TCGA
tcgafreq <- read.csv ( file = 'SOURCE_tcgafreqscreen_ALLTRANSCRIPTS.csv')
cclefreq <- read.csv ( file = 'SOURCE_cclefreqscreen_ALLTRANSCRIPTS.csv')
totalfreq = merge(tcgafreq, cclefreq, by = 'Var1', nomatch = 0)
write.csv (totalfreq, file ='outputs/totalfreq.csv')

# TRUNCATING PART
mutDATA <- read.csv (file = 'SOURCE_alltcga_jan2016.csv')
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


