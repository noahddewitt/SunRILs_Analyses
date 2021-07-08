library(LaByRInth)

setwd("/home/gbg_lab_admin/Array_60TB/Wheat_GBS/SunRILs_Mar2021")

filtVCFs <- list.files("populations")
filtVCFs <- filtVCFs[grepl("\\.vcf\\.gz$", filtVCFs)]

popInfo <- read.csv("xm_pedigrees.csv")

filtVCFs <- c("SunRILs_2021_UX2028_RILs_filt.vcf.gz", "SunRILs_2021_UX2029_RILs_filt.vcf.gz", "SunRILs_2021_UX2031_RILs_filt.vcf.gz")

for (vcf in filtVCFs) {
	print(paste0("Imputing file ", vcf, "..."))
	
	popName <- regmatches(vcf, regexpr("UX[0-9]{4}", vcf))
	parentNames <- as.character(popInfo[popInfo$Family == popName, c(2, 3)])

	vcf <- paste0("populations/", vcf)

	parentalResults <- LabyrinthImputeParents(
		vcf = vcf,
		out.file = paste0(popName, "_parents.rds"),
		parents = parentNames,
		generation = 5,
		geno.err = .02,
		parent.het = .001,
		parallel = TRUE,
		cores = 10)

	print(paste0("Parents ", paste0(parentNames, collapse = " and "), " imputed."))
	print("Imputing RILs...")

	popResult <- LabyrinthImputeProgeny(
		parental = parentalResults,
		out.file = paste0("populations/SunRILs_2021_", popName, "_RILs_imputed.vcf.gz"),
		parallel = TRUE,
		cores = 10)

	print(paste0("Population ", popName, " imputed."))
}
