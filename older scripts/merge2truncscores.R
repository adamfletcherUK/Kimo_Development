# merge two truncating scores

ccletrunc <- read.csv (file = 'SOURCE_lengthcorrectedfreq_ccleALL.csv')
tcgatrunc <- read.csv (file = 'SOURCE_lengthcorrectedfreq_tcgaALL.csv')

ccletruncscore <- as.data.frame(ccletrunc$Gene)
ccletruncscore$Score <- ccletrunc$score
colnames(ccletruncscore)[1] <- "Gene"

tcgatruncscore <- as.data.frame(tcgatrunc$Gene)
tcgatruncscore$Score <- tcgatrunc$score
colnames(tcgatruncscore)[1] <- "Gene"

ccletruncgenes <- ccletruncscore$Gene
tcgatruncgenes <- tcgatruncscore$Gene

unmatchedtruncsINTCGA <- setdiff (tcgatruncgenes, ccletruncgenes)
unmatchedtruncsINCCLE <- setdiff (ccletruncgenes, tcgatruncgenes)

mergedtruncs <- merge(tcgatruncscore, ccletruncscore, by = 'Gene',)
mergedtruncs$Score <- mergedtruncs$Score.x + mergedtruncs$Score.y
mergedtruncsfinal <-as.data.frame(mergedtruncs$Gene)
mergedtruncsfinal$Score <- mergedtruncs$Score
colnames(mergedtruncsfinal)[1] <- "Gene"

unmatchedtruncsINTCGA <- as.data.frame(unmatchedtruncsINTCGA)
colnames(unmatchedtruncsINTCGA)[1] <- "Gene"
unmatchedtruncsINTCGAmerge <-merge(unmatchedtruncsINTCGA, tcgatruncscore, by = 'Gene',)

unmatchedtruncsINCCLE <- as.data.frame(unmatchedtruncsINCCLE)
colnames(unmatchedtruncsINCCLE)[1] <- "Gene"
unmatchedtruncsINCCLEmerge <-merge(unmatchedtruncsINCCLE, ccletruncscore, by = 'Gene',)

totalScore <- rbind (mergedtruncsfinal, unmatchedtruncsINTCGAmerge, unmatchedtruncsINCCLEmerge)
write.csv (totalScore, file ='TEMP_mergedScores_largertranscript.csv')
