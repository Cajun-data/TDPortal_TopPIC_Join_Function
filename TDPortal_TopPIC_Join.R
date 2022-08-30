#' Read in TopPIC and TDportal files
#'
#'
#' @param outputname A character string specifying the file name to be exported
#'
#' @md
#'
#' @author James M Fulcher
#'
#' @export
#'
#' @import tidyverse
#' @import readxl
#' @importFrom dplyr "%>%"
#' 


library(tidyverse)
library(readxl)

join_results <- function(outputname)
  
{
  
tmp=list.files()
filenames1 <- tmp[grep(".xlsx", tmp)]

TDPortal <- lapply(filenames1, function(x) {readxl::read_excel(x)}) %>%
  do.call("rbind", .)

filenames2 <- tmp[grep(".RData", tmp)]

TopPIC <- lapply(filenames2, function(x){get(load(x,.GlobalEnv))}) %>%
  do.call("rbind", .)



### Change headers to match TDPortal
### Remove redundant or extra columns 
TopPIC <- TopPIC %>%
  dplyr::mutate(`Scan#` = `Scan(s)`) %>%
  dplyr::select(-`Scan(s)`) %>%
  dplyr::select(-ProteoForm, -IDs, -Condition, -AccMap, -AnnType,
         -`Special amino acids`, -CV, -`Retention time`, -Charge, -`Protein accession`,
         -`#matched peaks`, -`#matched fragment ions`, -isDecoy,-`#peaks`, -Condition,
         -`Data file name`)

datasets <- unique(TopPIC$Dataset)


TDPortal <- TDPortal %>%
  dplyr::mutate(Dataset = `File Name`,
                Dataset = gsub(".raw", "", Dataset),
                Dataset = substring(Dataset, 9)) %>%
  tidyr::separate_rows(`Fragment Scans`) %>%
  dplyr::mutate(`Scan#` = as.integer(`Fragment Scans`)) %>%
  dplyr::group_by(Dataset, `Scan#`) %>%
  dplyr::slice_min(order_by = `Global Q-value`, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-PFR, -`Uniprot Id`,
         -`% Cleavages`, -`P-score`, -`C-score`, -`File Name`, -HitId,
         -`Result Set`,
         -`Precursor Scans`) %>%
  dplyr::filter(Dataset %in% datasets)

names(TopPIC) <- paste0("TopPIC.", names(TopPIC))

names(TDPortal) <- paste0("TDPortal.", names(TDPortal))

combined <- dplyr::full_join(TopPIC, TDPortal, by = c("TopPIC.Dataset" = "TDPortal.Dataset",
                                                  "TopPIC.Scan#" = "TDPortal.Scan#"))


combined <- combined %>%
  dplyr::mutate(SeqConflict = case_when(TopPIC.cleanSeq != TDPortal.Sequence ~ "Yes",
                                 TopPIC.cleanSeq == TDPortal.Sequence ~ "No",
                                 TRUE ~ "NA"))

combined <- combined %>%
  dplyr::mutate(Software = case_when(is.na(`TDPortal.Proteoform Level`) ~ "TopPIC Only",
                              is.na(`TopPIC.Proteoform Level`) ~ "TDPortal Only",
                              SeqConflict != "NA" ~ "Both"))


write.csv(combined, file = outputname)

}

