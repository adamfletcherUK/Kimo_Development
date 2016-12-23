#ALT ALIGN SCRIPT - calculates conservation score as number of most frequent amino acid and mutation score is per kinase - need to edit table to only include regions with conservation score above 19
STK11 <- read.csv ( file = 'top30/STK11.csv')
MAP2K4 <- read.csv ( file = 'top30/MAP2K4.csv')
BMPR2 <- read.csv ( file = 'top30/BMPR2.csv')
ACVR2A <- read.csv ( file = 'top30/ACVR2A.csv')
TSSK1B <- read.csv ( file = 'top30/TSSK1B.csv')
STK32A <- read.csv ( file = 'top30/STK32A.csv')
STK24 <- read.csv ( file = 'top30/STK24.csv')
CSNK1A1 <- read.csv ( file = 'top30/CSNK1A1.csv')
MAP3K19 <- read.csv ( file = 'top30/MAP3K19.csv')
NEK5 <- read.csv ( file = 'top30/NEK5.csv')
CSNK1G3 <- read.csv ( file = 'top30/CSNK1G3.csv')
CHEK2 <- read.csv ( file = 'top30/CHEK2.csv')
MAP3K1 <- read.csv ( file = 'top30/MAP3K1.csv')
MAP3K12 <- read.csv ( file = 'top30/MAP3K12.csv')
RIPK3 <- read.csv ( file = 'top30/RIPK3.csv')
DYRK3 <- read.csv ( file = 'top30/DYRK3.csv')
NIM1K <- read.csv ( file = 'top30/NIM1K.csv')
LATS1 <- read.csv ( file = 'top30/LATS1.csv')
FER <- read.csv ( file = 'top30/FER.csv')
EPHA5 <- read.csv ( file = 'top30/EPHA5.csv')
CSNK2A1 <- read.csv ( file = 'top30/CSNK2A1.csv')
TGFBR2 <- read.csv ( file = 'top30/TGFBR2.csv')
NEK11 <- read.csv ( file = 'top30/NEK11.csv')
BRSK1 <- read.csv ( file = 'top30/BRSK1.csv')
EPHA7 <- read.csv ( file = 'top30/EPHA7.csv')
PRKG2 <- read.csv ( file = 'top30/PRKG2.csv')
PAK7 <- read.csv ( file = 'top30/PAK7.csv')
AMHR2 <- read.csv ( file = 'top30/AMHR2.csv')
GSK3B <- read.csv ( file = 'top30/GSK3B.csv')
CDKL3 <- read.csv ( file = 'top30/CDKL3.csv')


TOP30 <- rbind (STK11, MAP2K4, BMPR2, ACVR2A, TSSK1B, STK32A, STK24, CSNK1A1, MAP3K19, NEK5, CSNK1G3, CHEK2, MAP3K1, MAP3K12, RIPK3, DYRK3, NIM1K, LATS1, FER, EPHA5, CSNK2A1, TGFBR2, NEK11, BRSK1, EPHA7, PRKG2, PAK7, AMHR2, GSK3B, CDKL3)
write.csv (TOP30, file='TEMP_TOP30.csv')

#merge with genbank to get location of start of Catalytic frag
TOP30 <- merge(TOP30, genbankkinases, by = 'Gene')

#work out protein location of frag location
TOP30$Codon <- TOP30$Old_Loc + TOP30$loc_actual_gxg - 1

#MERGE WITH CCLE MUTATION DATA
CCLE    <- read.csv ( file = 'SOURCE_ccletrimmed_namechanged_HRT18removed_HEC1Bremoved.csv')
CCLE = subset(CCLE, mutation_type == 'Missense_Mutation')
CCLE$Ref <- str_extract(CCLE$amino_acid_change, "[A-Z]")
CCLE$Codon <- str_extract(CCLE$amino_acid_change, "[0-9]+")

mergedccle <- merge(CCLE, TOP30, by = c('Gene', 'Codon', 'Ref') )
mergedccle$pastedref <- paste(mergedccle$Gene, mergedccle$Codon, mergedccle$Tumor_Sample_Barcode, sep = "_" )
mergedccle$DUP= !duplicated(mergedccle$pastedref)
mergedccle <- subset(mergedccle, DUP == 'TRUE')

write.csv (mergedccle, file='outputs/TEMP_merged_align_ccle.csv')

#look for number of kinases altered per codon
mergedccle$pastedref2 <- paste(mergedccle$Align_Loc, mergedccle$Gene, sep = "_" )
mergedccle$DUP= !duplicated(mergedccle$pastedref2)
mergedcclepk <- subset(mergedccle, DUP == 'TRUE')

write.csv (mergedcclepk, file='outputs/TEMP_merged_align_perkinase_ccle.csv')

mergedcclepkTRIM <- as.data.frame (mergedcclepk$Gene)
mergedcclepkTRIM$Align_Loc <- mergedcclepk$Align_Loc
colnames(mergedcclepkTRIM)[1] <- "Gene"


#MERGE WITH TCGA MUTATION DATA
TCGA    <- read.csv ( file = 'SOURCE_allTCGA_jan2016.csv')
TCGA = subset(TCGA, mutation_type == 'Missense_Mutation')
TCGA$Ref <- str_extract(TCGA$amino_acid_change, "[A-Z]")
TCGA$Codon <- str_extract(TCGA$amino_acid_change, "[0-9]+")

mergedtcga <- merge(TCGA, TOP30, by = c('Gene', 'Codon', 'Ref') )
mergedtcga$pastedref <- paste(mergedtcga$Gene, mergedtcga$Codon, mergedtcga$case_id, sep = "_" )
mergedtcga$DUP= !duplicated(mergedtcga$pastedref)
mergedtcga <- subset(mergedtcga, DUP == 'TRUE')

write.csv (mergedtcga, file='outputs/TEMP_merged_align_tcga.csv')

#look for number of kinases altered per codon
mergedtcga$pastedref2 <- paste(mergedtcga$Align_Loc, mergedtcga$Gene, sep = "_" )
mergedtcga$DUP= !duplicated(mergedtcga$pastedref2)
mergedtcgapk <- subset(mergedtcga, DUP == 'TRUE')

write.csv (mergedtcgapk, file='outputs/TEMP_merged_align_perkinase_tcga.csv')

mergedtcgapkTRIM <- as.data.frame (mergedtcgapk$Gene)
mergedtcgapkTRIM$Align_Loc <- mergedtcgapk$Align_Loc
colnames(mergedtcgapkTRIM)[1] <- "Gene"



mergedBOTH <- rbind (mergedcclepkTRIM, mergedtcgapkTRIM)
mergedBOTH$pasted <- paste (mergedBOTH$Gene, mergedBOTH$Align_Loc, sep = "_")
mergedBOTH$DUP= !duplicated(mergedBOTH$pasted)
mergedBOTH <- subset(mergedBOTH, DUP == 'TRUE')


freq <-table(mergedBOTH$Align_Loc)
write.csv (freq, file ='outputs/TEMP_freq_perkinase_BOTH.csv')


#alt conservation score - calculate the score based on most common residue

conserve <- read.csv ( file = 'SOURCE_Top30_rotate.csv')

# first merge each position
conserve$merge <- paste (conserve$STK11, conserve$MAP2K4, conserve$BMPR2, conserve$ACVR2A, conserve$TSSK1B, conserve$STK32A, conserve$STK24, conserve$CSNK1A1, conserve$MAP3K19, conserve$NEK5, conserve$CSNK1G3, conserve$CHEK2, conserve$MAP3K1, conserve$MAP3K12, conserve$RIPK3, conserve$DYRK3, conserve$NIM1K, conserve$LATS1, conserve$FER, conserve$EPHA5, conserve$CSNK2A1, conserve$TGFBR2, conserve$NEK11, conserve$BRSK1, conserve$EPHA7, conserve$PRKG2, conserve$PAK7, conserve$AMHR2, conserve$GSK3B, conserve$CDKL3, sep = ",")
conserve$dash <- str_count(conserve$merge, pattern = "-")
conserve$A <- str_count(conserve$merge, pattern = "A")
conserve$C <- str_count(conserve$merge, pattern = "C")
conserve$D <- str_count(conserve$merge, pattern = "D")
conserve$E <- str_count(conserve$merge, pattern = "E")
conserve$F <- str_count(conserve$merge, pattern = "F")
conserve$G <- str_count(conserve$merge, pattern = "G")
conserve$H <- str_count(conserve$merge, pattern = "H")
conserve$I <- str_count(conserve$merge, pattern = "I")
conserve$K <- str_count(conserve$merge, pattern = "K")
conserve$L <- str_count(conserve$merge, pattern = "L")
conserve$M <- str_count(conserve$merge, pattern = "M")
conserve$N <- str_count(conserve$merge, pattern = "N")
conserve$P <- str_count(conserve$merge, pattern = "P")
conserve$Q <- str_count(conserve$merge, pattern = "Q")
conserve$R <- str_count(conserve$merge, pattern = "R")
conserve$S <- str_count(conserve$merge, pattern = "S")
conserve$T <- str_count(conserve$merge, pattern = "T")
conserve$V <- str_count(conserve$merge, pattern = "V")
conserve$W <- str_count(conserve$merge, pattern = "W")
conserve$Y <- str_count(conserve$merge, pattern = "Y")


conserve$score <- conserve$A
conserve$score <- ifelse (conserve$C < conserve$score, conserve$score, conserve$C)
conserve$score <- ifelse (conserve$D < conserve$score, conserve$score, conserve$D)
conserve$score <- ifelse (conserve$E < conserve$score, conserve$score, conserve$E)
conserve$score <- ifelse (conserve$F < conserve$score, conserve$score, conserve$F)
conserve$score <- ifelse (conserve$G < conserve$score, conserve$score, conserve$G)
conserve$score <- ifelse (conserve$H < conserve$score, conserve$score, conserve$H)
conserve$score <- ifelse (conserve$I < conserve$score, conserve$score, conserve$I)
conserve$score <- ifelse (conserve$K < conserve$score, conserve$score, conserve$K)
conserve$score <- ifelse (conserve$L < conserve$score, conserve$score, conserve$L)
conserve$score <- ifelse (conserve$M < conserve$score, conserve$score, conserve$M)
conserve$score <- ifelse (conserve$N < conserve$score, conserve$score, conserve$N)
conserve$score <- ifelse (conserve$P < conserve$score, conserve$score, conserve$P)
conserve$score <- ifelse (conserve$Q < conserve$score, conserve$score, conserve$Q)
conserve$score <- ifelse (conserve$R < conserve$score, conserve$score, conserve$R)
conserve$score <- ifelse (conserve$S < conserve$score, conserve$score, conserve$S)
conserve$score <- ifelse (conserve$T < conserve$score, conserve$score, conserve$T)
conserve$score <- ifelse (conserve$V < conserve$score, conserve$score, conserve$V)
conserve$score <- ifelse (conserve$W < conserve$score, conserve$score, conserve$W)
conserve$score <- ifelse (conserve$Y < conserve$score, conserve$score, conserve$Y)

write.csv (conserve, file='outputs/conserved.csv')

#combine scores

freqBOTH <- read.csv ( file = 'outputs/TEMP_freq_perkinase_BOTH.csv')
colnames(freqBOTH)[2] <- "Align_Loc"

comboscore <- merge(freqBOTH, conserve, by = 'Align_Loc' )
write.csv (comboscore, file='outputs/TEMP_comboscore.csv')
