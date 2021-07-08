#!/bin/bash

popNames=`cat sunRILs_popNames.txt`
PL_DIR="/home/gbg_lab_admin/Array_60TB/git_repos/TASSEL_GBS_pipeline"

parentKeyFile=`cat SunRILs_keyFile.txt | grep -Ev "UX[0-9]*-[0-9]*|RIL" | grep -v "species\|blank\|Blank\|BLANK\|suspect"`


for pop in $popNames
	do echo "Starting SNP calling for pop ${pop}" 
	cd "/home/gbg_lab_admin/Array_60TB/Wheat_GBS/SunRILs_Mar2021"

	#We also want the names of specifically the two parents used to generate the pop
	parentOne=`cat xm_pedigrees.csv | grep "${pop}" | cut -d, -f2`
	parentTwo=`cat xm_pedigrees.csv | grep "${pop}" | cut -d, -f3`

	
	mkdir "${pop}_productionDiscovery"
	WKDIR="/home/gbg_lab_admin/Array_60TB/Wheat_GBS/SunRILs_Mar2021/${pop}_productionDiscovery"

	#This will manualy keep these two parents in filtering step
	echo "${parentOne}" > "${WKDIR}/popParents.txt"
	echo "${parentTwo}" >> "${WKDIR}/popParents.txt"


	cat SunRILs_keyFile.txt | grep "species\|$pop" | grep -v "blank\|Blank\|BLANK" \
			       > "${pop}_productionDiscovery/SunRILs_${pop}_keyFile.txt"
	echo "$parentKeyFile" >> "${pop}_productionDiscovery/SunRILs_${pop}_keyFile.txt"

	cd "$PL_DIR"

	#Run discovery/production on parents + population
	bash tassel_disc_plus_prod_bwa.sh --workdir="${WKDIR}" \
                                       --study="SunRILs_2021_${pop}" \
                                       --keyfile="SunRILs_${pop}_keyFile.txt" \
                                       --fastq=/home/gbg_lab_admin/Array_60TB/fastq/Wheat/ \
                                       --ref=/home/gbg_lab_admin/Array_60TB/GBS_Reference_Genomes/Ensembl_v41_IWGSC_v1.0/Triticum_aestivum.IWGSC.dna.toplevel.fa \
                                       --enzymes=PstI-MspI \
                                       --taglength=85 \
				       --run-prod

	#Split off parents and population files. Parent files will be combined to capture all SNP

	cd "${WKDIR}/raw_VCF"

	gunzip "SunRILs_2021_${pop}_production.vcf.gz"

	#Lay down header for two new vcf files
	cat "SunRILs_2021_${pop}_production.vcf" | grep "##" > "SunRILs_2021_${pop}_RILs.vcf"
	cp "SunRILs_2021_${pop}_RILs.vcf" "SunRILs_2021_${pop}_parents.vcf"

#	popList=`cat ${WKDIR}/SunRILs_UX1990_keyFile.txt | cut -d$'\t' -f12 | grep -E "^UX"`
	#First, we identify the row position of all the biparental lines
	popIndices=`cat summary4.txt | cut -d$'\t' -f2 | grep -En "^UX|${parentOne}|${parentTwo}" | cut -d: -f1`
        parIndices=`cat summary4.txt | cut -d$'\t' -f2 | grep -Env "^UX|Taxa" | cut -d: -f1`
	#Then, we correct for the number of static VCF columns minus the header row of the keyfile and add in commas 
	popIndices=`echo $popIndices | sed -re 's/[0-9]+/$((\0+8))/g; s/^/echo /e; s/ /,/g'`
        parIndices=`echo $parIndices | sed -re 's/[0-9]+/$((\0+8))/g; s/^/echo /e; s/ /,/g'`

	cat "SunRILs_2021_${pop}_production.vcf" | grep -v "##" | cut -d$'\t' -f"1,2,3,4,5,6,7,8,9,${popIndices}"  >> "SunRILs_2021_${pop}_RILs.vcf"
	cat "SunRILs_2021_${pop}_production.vcf" | grep -v "##" | cut -d$'\t' -f"1,2,3,4,5,6,7,8,9,${parIndices}" >> "SunRILs_2021_${pop}_parents.vcf"

	#gzip "SunRILs_2021_${pop}_production.vcf"
	mv "SunRILs_2021_${pop}_parents.vcf" "${WKDIR}/../parents/" 

	cd "$PL_DIR"
	bash vcf_filter_parallel.sh --workdir="${WKDIR}" \
                   --vcfin="raw_VCF/SunRILs_2021_${pop}_RILs.vcf" \
                   --vcfout="filt_and_imp_VCF/SunRILs_2021_${pop}_RILs_filt.vcf.gz" \
                   --taxamiss=0.85 \
                   --maxdep=100 \
                   --snpmiss=0.3 \
                   --snphet=0.2 \
                   --maf=0.3 \
                   --removechr=UN \
                   --ncores=8 

	mv "${WKDIR}/raw_VCF/SunRILs_2021_${pop}_RILs.vcf" "${WKDIR}/../populations/"
	mv "${WKDIR}/filt_and_imp_VCF/SunRILs_2021_${pop}_RILs_filt.vcf.gz" "${WKDIR}/../populations/"

	echo "$pop" >>  "${WKDIR}/../sunRILs_completedPops.txt"

	echo "Finished SNP calling and filtering"
done

