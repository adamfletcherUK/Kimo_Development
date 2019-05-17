#Alignment script
#create catalytic frags for alignment
top20 <- read.csv ('SOURCE_top20truncs.csv')
catFrag <- merge(genbankkinases, top20, by = 'Gene', nomatch = 0)
catFrag$endOfFrag <- catFrag$loc_actual_APE + 40
catFrag$catFrag <- str_sub(catFrag$Protein_Seq, catFrag$loc_actual_gxg, catFrag$endOfFrag)
write.csv (catFrag, file ='output/catFrag.csv')