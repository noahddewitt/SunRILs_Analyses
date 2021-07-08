setwd("~/Array_60TB/Wheat_GBS/SunRILs_Mar2021")

library(gaston)
library(parallel)

numCores = 8

parentGenos <- read.csv("comb_parents_geno.csv", row.names = 1)
colnames(parentGenos) <- gsub("\\.", "-", colnames(parentGenos))
#colnames(parentGenos) <- gsub("LA0964C-P2", "LA09264C-P2", colnames(parentGenos))
snpInfo <- parentGenos[,c(1:2)]
parentGenos <- parentGenos[,-c(1:2)]

famPeds <- read.csv("xm_pedigrees.csv", stringsAsFactors = F) 

fileList <- list.files("populations")
filtVCFs <- fileList[grepl("RILs_imputed.vcf.gz", fileList)]

#Takes a cross name and SNP name, then decide pop genotypes (all 0, all 2, all NA). 
returnParentGeno <- function(snpName, famName, pedDat, parGeno) {
	famParents <- famPeds[famPeds$Family == famName, c(2,3)]	

	parOneGeno <- parGeno[snpName, famParents[1, 1]]
	parTwoGeno <- parGeno[snpName, famParents[1, 2]]
	
	if (length(parOneGeno) == 0) {
		print("Null result in returnParentGeno")
		print(parOneGeno)
		print(snpName)
		print(famParents)
		break
	}
	
	if ((is.na(parOneGeno)) | (parOneGeno == "1")) {
		if (parTwoGeno %in% c("0", "2")) {
				return(parTwoGeno)
		} else {
				return(NA)
		}
	} else if ((is.na(parTwoGeno)) | (parTwoGeno) == "1") {
		if (parOneGeno %in% c("0", "2")) {
				return(parOneGeno)
		} else {
				return(NA)
		}

	} else if ((parOneGeno == parTwoGeno) & (parOneGeno %in% c("0", "2"))) {
		return(parOneGeno)
	} else {
		return(NA)	
	}

}

#"inflates" vector of parental SNP calls into matrix for all indvs
inflateSNP <- function(snpMat, famVec) {
	inflatedMat <- NULL
	for (i in c(1:nrow(snpMat))) {
		subFamVec <- famVec[famVec == snpMat[i, 1]]
		famMat <- matrix(rep(snpMat[i,], length(subFamVec)), 
			nrow = length(subFamVec), byrow = T)
		colnames(famMat) <- colnames(snpMat) 
		rownames(famMat) <- subFamVec
		inflatedMat <- rbind(inflatedMat, famMat)
	}
	return(inflatedMat)
}

combVCF <- NULL
for (file in filtVCFs) {
	print(paste0("Merging file ", file))
	curVCF <- read.vcf(paste0("populations/", file), convert.chr = F)
	curVCF <- as.matrix(curVCF)

	curVCF <- curVCF[grepl("^UX", rownames(curVCF)), ]

	#Modify file to contain cross info column. Will be used now, and also to later expand 
	#Larger data frame. Ie, both dataframes are expanded columnwise, then rbinded. 
	famVec <- gsub("-.*$", "", rownames(curVCF))
	#Sanity check
	if (length(unique(famVec)) == 1) {
		curVCF <- cbind(Family = famVec, curVCF)	
	} else {
		print("Multiple families in VCF file!")
		break
	}

	#In about 20% of SNP where SNP was segging in pop, both parents had same 
	#Homozygous call. But that's okay, because we're not imputing those anyways.
	#But it could be indicative of issues with samples -- thankfully we'ere getting 
	#More samples. 

	absentSNP <- colnames(combVCF)[!colnames(combVCF) %in% colnames(curVCF)]
	newSNP <- colnames(curVCF)[!colnames(curVCF) %in% colnames(combVCF)] #Fam?

	absentSNPCalls <- mclapply(absentSNP, returnParentGeno, 
			famName = famVec[1], pedDat = famPeds, parGeno = parentGenos,
			mc.cores = numCores)
	absentSNPCalls <- unlist(absentSNPCalls)
	names(absentSNPCalls) <- absentSNP
	absentSNPCalls <- as.data.frame(t(c(Family = famVec[1], absentSNPCalls)))

	absentSNPFrame <- inflateSNP(absentSNPCalls, curVCF[,1])
	curVCF <- cbind(curVCF, absentSNPFrame[,-1])

	if (!is.null(combVCF)) {
		newSNPCalls <- NULL
		for (curFam in unique(combVCF[,1])) {
			newSNPCallsFam <- mclapply(newSNP, returnParentGeno,
					famName = curFam, pedDat = famPeds, 
					parGeno = parentGenos, mc.cores = numCores)
			newSNPCallsFam <- unlist(newSNPCallsFam)
			names(newSNPCallsFam) <- newSNP
			newSNPCalls <- rbind(newSNPCalls, c(Family = curFam, newSNPCallsFam))

		}

		newSNPFrame <- inflateSNP(newSNPCalls, combVCF[,1])
		combVCF <- cbind(combVCF, newSNPFrame[,-1])
	}

	colOrder <- c("Family", colnames(curVCF)[-1][order(colnames(curVCF)[-1])])
	curVCF <- curVCF[, colOrder]
	combVCF <- combVCF[, colOrder]

	combVCF <- rbind(combVCF, curVCF)	

}

#combVCF <- combVCF[, !grepl("^SUN", colnames(combVCF))]

#Otherwise sorts e.g. 800 after 6000
snpChr <- gsub("_\\d+$", "", colnames(combVCF)[-1])
snpPos <- as.numeric(gsub("S[0-9][A,B,D]_", "", colnames(combVCF)[-1]))

trueSNPOrder <- colnames(combVCF)[-1][order(snpChr, snpPos)] 
combVCF <- combVCF[,trueSNPOrder]

#Now we re-add back in the parents. This serves as the "filtering" stage for the parents,
#Where the only parental SNPs that are preserved passed filtering in at least one subpop

subParentGenos <- t(parentGenos[trueSNPOrder,])
combVCF <- rbind(subParentGenos, combVCF)

#Replace major/minor genotype coding with nucleotide coding

nucCombVCF <- NULL
for (snp in colnames(combVCF)) {
	genoVec <- combVCF[,snp]
	
	genoVec[genoVec == "0"] <- paste0(snpInfo[snp,1], snpInfo[snp,1])		
	genoVec[genoVec == "1"] <- paste0(snpInfo[snp,1], snpInfo[snp,2])		
	genoVec[genoVec == "2"] <- paste0(snpInfo[snp,2], snpInfo[snp,2])		
	genoVec[is.na(genoVec)] <- "NN"

	nucCombVCF <- cbind(nucCombVCF, genoVec)
}
colnames(nucCombVCF) <- colnames(combVCF)

#Tassel plugin wants the file in hapmap format.
combHmp <- cbind("rs#" = trueSNPOrder, alleles = c(NA), 
		chrom = gsub("^S", "", gsub("_\\d+$", "", trueSNPOrder)),
		pos = as.numeric(gsub("S[0-9][A,B,D]_", "", trueSNPOrder)),
		strand = c(NA), "assembly#" = c(NA), center = c("ERSGGL"),
		protLSID = c(NA), assayLSID = c("GBS"), panelLSID = c("SunRILs"),
		QCode = c(NA))  

combHmp <- cbind(combHmp, t(nucCombVCF))

write.table(combHmp, "SunRILs_combGenos.hmp.txt", row.names = F, sep = "\t")
