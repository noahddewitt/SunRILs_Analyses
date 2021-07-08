rm(list = ls())
library(tidyverse)
library(DBI)

setwd("/home/gbg_lab_admin/Array_60TB/Wheat_GBS/SunRILs_Mar2021")

con <- dbConnect(RMariaDB::MariaDB(), group = "client", dbname = "genotype")
res <- dbSendQuery(con, "SELECT g.species, g.hi_seq_library_id, h.Flowcell, 
  h.Lane, g.Barcode, g.FullSampleName, g.enzyme, g.study_code, g.library_prep_id, 
  g.plate_id, g.well, g.original_sample_name, g.study_name, g.note
  FROM genotype.isolations AS g
  INNER JOIN genotype.hi_seq_library_id AS h ON g.hi_seq_library_id = h.hi_seq_library_id
  ORDER BY g.library_prep_id;")
key <- dbFetch(res)
dbClearResult(res)
dbDisconnect(con)

key <- filter(key, str_detect(study_name, "XM|LA95135 x AGS2000|MPV57 x LA95135|LA95135 x MPV57"))

#Get rid of library Ids because only some have them, also we want to combine same genotypes
key <- mutate(key, FullSampleName = if_else(str_detect(FullSampleName, "_\\d{6}$"), 
str_replace(FullSampleName, "_\\d{6}$", "") , FullSampleName)) %>%
	filter(!str_detect(FullSampleName, "BLANK|Blank|blank|Water|WATER|water|NO-DNA|SUSPECT")) %>%
	mutate(study_name = str_replace(study_name, ".*UX", "UX"),
	       study_name = str_replace(study_name, "\\-4$", ""),
	       study_name = str_replace(study_name, "LA95135 x AGS2000", "UX1444")) %>%
	mutate(FullSampleName = str_replace(FullSampleName, "RIL-", "UX1443-"))

#All-number FullSampleNames are from the LA population (with one exception).
key <- mutate(key, FullSampleName = if_else((study_name == "UX1444") & (!is.na(as.numeric(FullSampleName))), paste0("UX1444-", FullSampleName), FullSampleName)) %>%
	mutate(FullSampleName = if_else(FullSampleName == "225-NWG", "UX1444-225-NWG", FullSampleName)) %>%
	mutate(FullSampleName = if_else(str_detect(FullSampleName, "^LM-"), paste0(str_replace(FullSampleName, "^LM-", ""), "-RESEQ"), FullSampleName))

write_tsv(key, "SunRILs_keyFile.txt")

#Now, want to also turn the key file into a pedigree file for 
