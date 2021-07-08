setwd("~/Array_60TB/Wheat_GBS/SunRILs_Mar2021/parents")

library(gaston)

combSNP <- NULL
for (file in list.files()) {
	curVCF <- read.vcf(file, convert.chr = F)
	snpInfo <- curVCF@snps
	#curVCF <- select.snps(curVCF, maf > .05) don't wanna filter
	curVCF <- as.matrix(curVCF)
	curVCF <- t(curVCF)

	#We assume that all the genos should be the same. 
	curVCF <- curVCF[!rownames(curVCF) %in% rownames(combVCF), ]
	snpInfo <- snpInfo[!snpInfo$id %in% rownames(combVCF),]

	curVCF <- cbind(A1 = snpInfo$A1, A2 = snpInfo$A2, curVCF)

	combVCF <- rbind(combVCF, curVCF)
}

combVCF <- combVCF[order(rownames(combVCF)), ]
#combVCF <- combVCF[, !grepl("^SUN", colnames(combVCF))]

write.csv(combVCF, "../comb_parents_geno.csv")
