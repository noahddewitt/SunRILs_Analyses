#SunRILs data processing and analysis

Connected SunRILs population constructed from shared parents representing elite southeastern soft red winter wheat germplasm lines. The actual data for this project is contained in the ERSGGL database, but this repo should serve as documentation of the pipeline used to assemble genotype data and the analyses performed using that data. 

##SNP Calling and Imputation

Calling all lines together and using beagle imputation worked pretty well in the initial test run. Theoretically though, as we have a set of related inbred biparental populations, treating each population separately should give us better power to call SNP private to individual populations, filter SNP, and perform imputation. The following scripts were run in the listed order:

xm\_pedigrees.csv -- Manually prepared csv file containing information on what parents were used to generate which populations, coded by the standard ERSGGL UX\_\_\_ cross format. 

makeKeyFile.R -- Tidy R script which usesRMariaDB to pull out all the records of the lines in the XM population contained on the server's SQL database, perform cleaning of line names for consistent formatting, and then writing the key file out in a format usable by TASSEL for SNP calling.

callIndividualPops.sh -- Bash script that iterates through populations and performs SNP calling and filtering using TASSEL and standard ERSGGL genotyping bash scripts. The larger key file is subsetted into all parents and the lines in each population, and a combined discovery and production run is performed for each. A VCF file of the parental genotypes for all SNP discovered by TASSEL is spit out as a record used to compare TASSEL runs in different populations, and a VCF file consisting of the RILs in that population and their parents are subject to filitering based on biparental population-specific criteria before being written as a VCF. 

labyrinthImputePops.R -- Base R script that performs imputation on the RILs of each population using the LB-Impute viterbi HMM algorithm as implemented in the R package LaByRInth. 

mergeParentVCFs.R -- Base R script that merges all of the parental VCF files produced by callIndividualPops.sh. This file serves as a record of the parental genotypes discovered in ANY of the different TASSEL runs.

mergePopulationVCFs.R -- Base R script that merges the SNP calls of all the different populations together into a single hapmap file suitable for downstream analyses. Given that each population segregates for only a portion of the total markers in all of the populations, and markers fixed within each population are not discovered by TASSEL, the script uses the parental records produced by the mergeParentVCFs.R script to impute fixed genotypes within each population from the parental genotypes. 


