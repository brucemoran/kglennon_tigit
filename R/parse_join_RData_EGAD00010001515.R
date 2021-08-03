#! R

##NB data is under a Data Access Committee and so is not shared
##This file is to show how data was parsed and saved into the object used for the DGE analysis

##parse files for counts, join and save as RData
library(tidyverse)
library(biomaRt)

samples <- dir(pattern = "EGAF")
master_tb <- tibble()

for(x in 1:length(samples)){
  print(paste0("Working on: ", samples[x]))

  rcc <- paste0(samples[x], "/", dir(samples[x], pattern = ".RCC$"))
  rcc_smp <- gsub(".RCC", "", paste(c("NS", unlist(strsplit(rcc, "_")[[1]])[3:5]), collapse = "_"))
  rcc_tb <- readr::read_csv(rcc, skip = grep("Code_Summary", unlist(readr::read_csv(rcc)))[1]+4) %>%
            dplyr::rename(!!quo_name(rcc_smp) := "Count")

  if(x == 1){
    master_tb <- rcc_tb
  } else {
    master_tb <- dplyr::left_join(master_tb, rcc_tb)
  }
}

##remove NA
master_tb <- master_tb %>% dplyr::filter(!is.na(Accession))
master_tb$Accession_nv <- unlist(lapply(master_tb$Accession, function(f){
    strsplit(f, "\\.")[[1]][1]
  }))

master_tb <- dplyr::select(.data = master_tb, 1, 2, 3, Accession_nv, everything())

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
refseq_mb_tb <- tibble::as_tibble(biomaRt::getBM(attributes=c("refseq_mrna",
                                       "external_gene_name",
                                       "transcript_length"),
                          mart = mart)) %>%
                dplyr::filter(refseq_mrna %in% master_tb$Accession_nv)
refseq_mbl_tb <- refseq_mb_tb %>% dplyr::group_by(refseq_mrna) %>%
                 dplyr::summarize(mean_length = mean(transcript_length)) %>%
                 dplyr::ungroup()

counts.mat <- master_tb %>% dplyr::filter(Accession_nv %in% refseq_mbl_tb$refseq_mrna) %>%
              dplyr::arrange(Accession_nv) %>%
              as.data.frame() %>%
              column_to_rownames("Accession_nv") %>%
              dplyr::select(-1,-2,-3) %>%
              as.matrix()

##https://support.bioconductor.org/p/91218/
cml <- counts.mat / unlist(refseq_mbl_tb$mean_length)
tpm.mat <- t( t(cml) * 1e6 / colSums(cml) )

EGAD00010001515_tpm_tb <- dplyr::left_join(refseq_mb_tb[,c(1,2)], as_tibble(tpm.mat, rownames = "Accession_nv"), by = c("refseq_mrna" = "Accession_nv")) %>% dplyr::distinct()

EGAD00010001515_log2tpm_tb <- rapply(EGAD00010001515_tpm_tb, f = function(ff){log2(ff+0.00001)}, classes = c("numeric", "integer"), how = "replace")

##############
## Metadata ##
##############
url_metadata <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867418304458-mmc2.xlsx"
temp_metadata <- tempfile()
utils::download.file(url_metadata, temp_metadata)

##https://stackoverflow.com/questions/40857694/remove-everything-after-the-last-underscore-of-a-column-in-r
metadata <- tibble::as_tibble(readxl::read_xlsx(path = temp_metadata,
                                              sheet = "Cohort")) %>%
            dplyr::mutate("ns_sample_id" = unlist(lapply(sample_id, function(f){
              paste0("NS_", f)
            }))) %>%
            dplyr::mutate(site = sub("_[^_]+$", "", site_id))

smp_metadata <- unlist(lapply(grep("NS_", colnames(EGAD00010001515_tpm_tb), value = TRUE), function(f){
    smp <- paste(strsplit(f, "_")[[1]][2:3], collapse = "_")
    if(smp %in% metadata$sample_id){
      return(smp)
    }
  })) %>% sort() %>% unique()

EGAD00010001515_metadata <- metadata %>% dplyr::filter(sample_id %in% smp_metadata)

##save
EGAD00010001515_raw_count_tb <- master_tb
save(EGAD00010001515_raw_count_tb, EGAD00010001515_tpm_tb, EGAD00010001515_log2tpm_tb, EGAD00010001515_metadata, file = "EGAD00010001515.RData")
