#!/bin/bash

popNames=`cat ../sunRILs_popNames.txt`
PL_DIR="/home/gbg_lab_admin/Array_60TB/git_repos/TASSEL_GBS_pipeline"

for pop in $popNames
	do echo "Starting Sungrains SNP calling based on pop ${pop} db" 

	mkdir "${pop}_sungrains_production"
	WKDIR="/home/gbg_lab_admin/Array_60TB/Wheat_GBS/SunRILs_Mar2021/diverseproductions/${pop}_sungrains_production"

	cd "$PL_DIR"
	bash tassel_production.sh --workdir="${WKDIR}" \
                               --database=/home/gbg_lab_admin/Array_60TB/Wheat_GBS/SunRILs_Mar2021/"${pop}"_productionDiscovery/database/SunRILs_2021_"${pop}".db \
                               --study=Sungrains_"{pop}" \
                               --keyfile=/home/gbg_lab_admin/Array_60TB/Wheat_GBS/SUNGRAINS_May2021/corrected_merged_key.txt \
                               --fastq=/home/gbg_lab_admin/Array_60TB/fastq/Wheat/ \
                               --ref=/home/gbg_lab_admin/Array_60TB/GBS_Reference_Genomes/Ensembl_v41_IWGSC_v1.0/Triticum_aestivum.IWGSC.dna.toplevel.fa \
                               --enzymes='PstI-MspI' \
                               --taglength=85


	## Filter the VCF file
	cd "$PL_DIR"
	bash vcf_filter_parallel.sh --workdir="${WKDIR}" \
                   --vcfin=raw_VCF/Sungrains_"${pop}"_production.vcf.gz \
                   --vcfout=filt_and_imp_VCF/Sungrains_"${pop}"_preimp_filt.vcf.gz \
                   --taxamiss=0.85 \
                   --taxahet=0.3 \
                   --maxdep=100 \
                   --snpmiss=0.2 \
                   --snphet=0.1 \
                   --maf=0.05 \
                   --removechr=UN \
                   --ncores=7

	## Imputation with Beagle 5.1
	cd "$PL_DIR"/imputation
	bash simple_beagle_imputation_w_map.sh "${WKDIR}"/filt_and_imp_VCF/Sungrains_"{pop}"_preimp_filt.vcf.gz SynOp_RIL906_v1.0_GBS_monotonic.map 6

	## Post-imputation filtering
	cd ..
	bash vcf_filter_parallel.sh --workdir="${WKDIR}" \
                   --vcfin=filt_and_imp_VCF/Sungrains_"${pop}"_preimp_filt_imp.vcf.gz \
                   --vcfout=postimp_filt_VCF/Sungrains_"${pop}"_postimp_filt.vcf.gz \
                   --maf=0.05 \
                   --snphet=0.1 \
                   --ncores=7


	echo "$pop" >>  "${WKDIR}/../Sungrains_completedPops.txt"

	echo "Finished SNP calling and filtering for pop ${pop}"
done

