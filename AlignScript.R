#NEW ALIGN SCRIPT
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

freq <-table(mergedcclepk$Align_Loc)
write.csv (freq, file ='outputs/TEMP_freq_perkinase_ccle.csv')

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

freq <-table(mergedtcgapk$Align_Loc)
write.csv (freq, file ='outputs/TEMP_freq_perkinase_tcga.csv')


# to make a freq table for both ccle and tcga need to merge prior to making a freq table

mergedccletrim <- as.data.frame (mergedcclepk$Gene)
mergedccletrim$Align_Loc <- mergedcclepk$Align_Loc
colnames(mergedccletrim)[1] <- "Gene"

mergedtcgatrim <- as.data.frame (mergedtcgapk$Gene)
mergedtcgatrim$Align_Loc <- mergedtcgapk$Align_Loc
colnames(mergedtcgatrim)[1] <- "Gene"

mergedboth <- rbind (mergedccletrim, mergedtcgatrim)
mergedboth$pastedref3 <- paste(mergedboth$Align_Loc, mergedboth$Gene, sep = "_" )
mergedboth$DUP= !duplicated(mergedboth$pastedref3)
mergedboth <- subset(mergedboth, DUP == 'TRUE')

write.csv (mergedboth, file='outputs/TEMP_merged_align_perkinase_both.csv')

freq <-table(mergedboth$Align_Loc)
write.csv (freq, file ='outputs/TEMP_freq_perkinase_both.csv')


#now work out a conservation score

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


conserve$A <- ifelse(conserve$A <1, 0, 1)
conserve$C <- ifelse(conserve$C <1, 0, 1)
conserve$D <- ifelse(conserve$D <1, 0, 1)
conserve$E <- ifelse(conserve$E <1, 0, 1)
conserve$F <- ifelse(conserve$F <1, 0, 1)
conserve$G <- ifelse(conserve$G <1, 0, 1)
conserve$H <- ifelse(conserve$H <1, 0, 1)
conserve$I <- ifelse(conserve$I <1, 0, 1)
conserve$K <- ifelse(conserve$K <1, 0, 1)
conserve$L <- ifelse(conserve$L <1, 0, 1)
conserve$M <- ifelse(conserve$M <1, 0, 1)
conserve$N <- ifelse(conserve$N <1, 0, 1)
conserve$P <- ifelse(conserve$P <1, 0, 1)
conserve$Q <- ifelse(conserve$Q <1, 0, 1)
conserve$R <- ifelse(conserve$R <1, 0, 1)
conserve$S <- ifelse(conserve$S <1, 0, 1)
conserve$T <- ifelse(conserve$T <1, 0, 1)
conserve$V <- ifelse(conserve$V <1, 0, 1)
conserve$W <- ifelse(conserve$W <1, 0, 1)
conserve$Y <- ifelse(conserve$Y <1, 0, 1)

conserve$total <- as.numeric(conserve$A + conserve$C + conserve$D + conserve$E + conserve$F + conserve$G + conserve$H + conserve$I + conserve$K + conserve$L + conserve$M + conserve$N + conserve$P + conserve$Q + conserve$R + conserve$S + conserve$T + conserve$V + conserve$W + conserve$Y)

conserve$score <- ifelse(conserve$dash >5, 20, conserve$total)

conserve$score <- 20 - conserve$score

write.csv (conserve, file='outputs/conserved.csv')

#merge conserved and freq (rougher version)
freqboth <- read.csv ( file = 'TEMP_freq.csv')
cons <- read.csv ( file = 'TEMP_conserved_alt.csv')
mergedscore <- merge(freqboth, cons, by = 'Align_Loc')
write.csv (mergedscore, file='outputs/mergedscore.csv')

# merge with frequency
#conserve <- read.csv ( file = 'SOURCE_conserved_scores.csv')
#align_freq <- read.csv ( file = 'SOURCE_freq_align_mutations_both.csv')
#combinedscore <- merge(conserve, align_freq, by = 'Align_Loc' )
#combinedscore$combo <- as.numeric (combinedscore$Freq * combinedscore$score)
#write.csv (combinedscore, file='output/combinedscore.csv')

#more stringent version
#conserve_merge <- conserve$merge
#conserve <- as.data.frame(conserve$Align_Loc)
#conserve$Merge <- conserve_merge
#colnames(conserve)[1] <- "Align_Loc"

stringentTCGA <- merge (mergedtcga, conserve, by = 'Align_Loc')
write.csv (stringentTCGA, file ='TEMP_check.csv')
stringentTCGA <- stringentTCGA[,c("Align_Loc","Gene","Ref", "case_id", "amino_acid_change", "merge")]

stringentCCLE <- merge (mergedccle, conserve, by = 'Align_Loc')
stringentCCLE <- stringentCCLE[,c("Align_Loc","Gene","Ref", "Tumor_Sample_Barcode", "amino_acid_change", "merge")]
colnames(stringentCCLE)[4] <- "case_id"
write.csv (stringentCCLE, file ='TEMP_check.csv')

stringentCOMBO <- rbind (stringentCCLE, stringentTCGA)

write.csv (stringentCOMBO, file ='TEMP_check.csv')
#identify VAR
stringentCOMBO$Var <- str_sub(stringentCOMBO$amino_acid_change, -1, -1 )

#locate matches
stringentCOMBO$Present  <- str_count(stringentCOMBO$merge,  stringentCOMBO$Var)

write.csv (stringentCOMBO , file='outputs/STRINGENTcombinedscore.csv')
# subset only mutations with Present score of 0 (i.e. var not seen naturally in 20 kinases)
stringentMUTATIONS <- subset (stringentCOMBO, Present == '0')
#remove gene duplications between ccle and tcga
stringentMUTATIONS$pastedref <- paste(stringentMUTATIONS$Align_Loc, stringentMUTATIONS$Gene, sep = "_" )
stringentMUTATIONS$DUP= !duplicated(stringentMUTATIONS$pastedref)
stringentMUTATIONS <- subset(stringentMUTATIONS, DUP == 'TRUE')
write.csv (stringentMUTATIONS , file='outputs/stringentMUTATIONS.csv')
#frequency
freq <-table(stringentMUTATIONS$Align_Loc)
write.csv (freq, file ='TEMP_freq_STRINGENT.csv')

#merge conserved and freq (rougher version)
freqboth <- read.csv ( file = 'TEMP_freq_STRINGENT.csv')
cons <- read.csv ( file = 'TEMP_conserved_alt.csv')
mergedscore <- merge(freqboth, cons, by = 'Align_Loc')
write.csv (mergedscore, file='outputs/mergedscore.csv')


# merge with stringent with conservation
conserve <- read.csv ( file = 'SOURCE_conserved_scores.csv')
stringent_freq <- read.csv ( file = 'SOURCE_freq_stringent_mutations.csv')
combinedSTRINGENTscore <- merge(conserve, stringent_freq, by = 'Align_Loc' )
combinedSTRINGENTscore$combo <- as.numeric (combinedSTRINGENTscore$Freq * combinedSTRINGENTscore$score)
write.csv (combinedSTRINGENTscore, file='output/STRINGENTcombinedscore.csv')

