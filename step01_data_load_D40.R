# 
# step01_data_load.R
# 

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("/home/swu/biostat_sns_20230207_rpt006876_g360_vcv2_performance_control_verification_report-main/")

# load libraries
library(tidyverse)
library(readr)
library(pbapply)
library(openxlsx)
library(RCurl)
library(jsonlite)

# source column types for ease of data read in
source("Code/col_types.R")

# read in negative site data
snv_indel_cnv_negatives <- openxlsx::read.xlsx("RawData/vcv2_negative_sites.xlsx", sheet = "snv_indel_cnv") %>%
  mutate(chrom = as.character(chrom))
fusion_negatives <- openxlsx::read.xlsx("RawData/vcv2_negative_sites.xlsx", sheet = "fusion")

# get batch ID to retrieve FCID and data for analysis
#batchid <- read_csv("RawData/batchid_g360_vcv2_control_performance_verification.csv", 
#                    col_types = col_types_batch)

# locate FC with wtf tool
#fc_dir <- RCurl::getURL(paste0("https://wtf.ghdna.io/flowcell/output/", batchid$FlowcellID), 
#                        ssl.verifyhost=FALSE,
#                        ssl.verifypeer=FALSE)

#if input FCID manually by SFW
#fc_dir <- c("/ghds/ivd/flowcentral/230130_NB551031_0031_AHTNYHBGXN.rundeck2",
#            "/ghds/ivd/flowcentral/230130_NB552393_0080_AHTKGGBGXN.rundeck1",
#            "/ghds/ivd/flowcentral/230203_NB552393_0081_AHTTTWBGXN.rundeck1",
#            "/ghds/ivd/flowcentral/230203_NB552398_0092_AHTKMTBGXN.rundeck1",
#            "/ghds/ivd/flowcentral/230211_NB552398_0093_AHVKNKBGXN.rundeck2",
#            "/ghds/ivd/flowcentral/230211_NB551031_0034_AHVJMKBGXN.rundeck3")
#fc_dir <- jsonlite::fromJSON(fc_dir)$sentinel[[1]]$path

#fc_dir1 <- "/ghds/ivd/flowcentral/240112_NB551031_0049_AHTNFCBGXN"
fc_dir2 <- "/ghds/ivd/flowcentral/240123_NB551031_0050_AHTMKKBGXN"


# get sample ids

fc_dir <- c(fc_dir2)
variant_agreement <- map(fc_dir, ~{
  sample_ids <- read_tsv(file.path(.x, "autoqc_sample_qc.hdr.tsv"), 
                         col_types = col_types_autoqc) %>%
    filter(category == "sample",
           grepl("VCv2", run_sample_id)) %>%
    pull(run_sample_id) %>%
    unique()
  
  variant_data <- pblapply(sample_ids, function(i){
    
    # i_snv     <- read_tsv(file.path(.x, paste0(i, ".snv_call.hdr.tsv")), col_types = col_types_snv)
    # i_indel   <- read_tsv(file.path(.x, paste0(i, ".indel_call.hdr.tsv")), col_types = col_types_indel)
    # i_cnv     <- read_tsv(file.path(.x, paste0(i, ".cnv_call.hdr.tsv")), col_types = col_types_cnv)
    # i_fusion  <- read_tsv(file.path(.x, paste0(i, ".fusion_call.hdr.tsv")), col_types = col_types_fusion)
    # i_msi  <- read_tsv(file.path(.x, paste0(i, ".msi_parameters.hdr.tsv")), col_types = col_types_msi)
    # i_ghm <- read_tsv(file.path(.x, paste0(i, ".ghm_snp.hdr.tsv")), col_types = col_types_ghm_snp)
    i_snv     <- data.table::fread(file.path(.x, paste0(i, ".snv_call.hdr.tsv")), colClass = c("chrom" = "character",
                                                                                               "run_sample_id" = "character",
                                                                                               "gene" = "character",
                                                                                               "mut_nt" = "character"))
    i_indel   <- data.table::fread(file.path(.x, paste0(i, ".indel_call.hdr.tsv")), colClass = c("chrom" = "character",
                                                                                                 "run_sample_id" = "character",
                                                                                                 "gene" = "character",
                                                                                                 "mut_nt" = "character"))
    i_cnv     <- data.table::fread(file.path(.x, paste0(i, ".cnv_call.hdr.tsv")), colClass = c("chrom" = "character",
                                                                                               "run_sample_id" = "character",
                                                                                               "gene" = "character"))
    i_fusion  <- data.table::fread(file.path(.x, paste0(i, ".fusion_call.hdr.tsv")), colClass = c("run_sample_id" = "character"))
    i_msi  <- data.table::fread(file.path(.x, paste0(i, ".msi_parameters.hdr.tsv")), colClass = c("chrom" = "character",
                                                                                                  "run_sample_id" = "character",
                                                                                                  "gene" = "character"))
    i_ghm <- data.table::fread(file.path(.x, paste0(i, ".ghm_snp.hdr.tsv")), colClass = c("chrom" = "character",
                                                                                          "run_sample_id" = "character"))
    
  
    # SNV positives
    i_snv_pos <- i_snv %>%
      filter((gene == "EGFR" & mut_aa == "L858R" & position == 55259515 & mut_nt == "T>G") | 
               (gene == "EGFR" & mut_aa == "T790M" & position == 55249071 & mut_nt == "C>T") | 
               (gene == "MET" & mut_aa == "" & position == 116412044 & mut_nt == "G>T")) %>%
      mutate(putative_call = "positive")
    
    # SNV negatives
    i_snv_neg <- snv_indel_cnv_negatives %>%
      filter(class == "snv") %>%
      mutate(run_sample_id = i) %>% 
      left_join(i_snv, by = c("run_sample_id", "gene", "chrom", "position", "mut_nt")) %>%
      mutate(putative_call = "negative")
    
    
    # indel
    i_indel_pos <- i_indel %>%
      filter((gene == "TP53" & mut_aa == "p.Leu35fs" & position == 7579582 & mut_nt == "C>CA") |  
               (gene == "EGFR" & mut_aa == "p.Glu746_Ala750del" & position == 55242465 & mut_nt == "GGAATTAAGAGAAGCA>G")) %>%
      mutate(putative_call = "positive")
    
    i_indel_neg <- snv_indel_cnv_negatives %>%
      filter(class == "indel") %>%
      mutate(run_sample_id = i) %>% 
      left_join(i_indel, by = c("run_sample_id", "gene", "chrom", "position", "mut_nt")) %>%
      mutate(putative_call = "negative")
    
    # CNV
    i_cnv_pos <- i_cnv %>%
      filter(gene %in% c("ERBB2", "MET")) %>%
      mutate(putative_call = "positive")
    
    i_cnv_neg <- snv_indel_cnv_negatives %>%
      filter(class == "cnv") %>%
      select(gene) %>%
      mutate(run_sample_id = i) %>% 
      left_join(i_cnv, by = c("run_sample_id", "gene")) %>%
      mutate(putative_call = "negative")
    
    # fusions
    i_fusion_pos  <- i_fusion %>%
      filter((gene_a == "ALK" & gene_b == "EML4" & downstream_gene == "A" & pos_a == 29448043 & pos_b == 42527768) | 
               (gene_a == "EML4" & gene_b == "ALK" & downstream_gene == "B" & pos_b == 29448043 & pos_a == 42527768)) %>%
      mutate(putative_call = "positive")
    
    i_fusion_neg <- fusion_negatives %>%
      select(gene_a, gene_b) %>%
      mutate(run_sample_id = i) %>% 
      left_join(i_fusion, by = c("run_sample_id", "gene_a", "gene_b")) %>%
      mutate(putative_call = "negative")
    
    # msi
    i_msi_pos <- i_msi %>%
      filter(gene == "ARID1A", grepl("REPEAT_SEQ=GCA", notes), grepl("LENGTH=8", notes), start_pos == 27100181, chrom == 1) %>%
      mutate(putative_call = "positive")
    i_msi_neg <- i_msi %>%
      filter((gene == "IPO5"  & grepl("REPEAT_SEQ=T", notes) & grepl("LENGTH=8", notes) & chrom == 13 & start_pos == 98668972) |
               (gene == "BRCA2" & grepl("REPEAT_SEQ=A", notes) & grepl("LENGTH=7", notes) & chrom == 13 & start_pos == 32953632) |
               (gene == "UBR5"  & grepl("REPEAT_SEQ=A", notes) & grepl("LENGTH=7", notes) & chrom == 8  & start_pos == 103289457) |
               (gene == "FGFR2" & grepl("REPEAT_SEQ=C", notes) & grepl("LENGTH=7", notes) & chrom == 10 & start_pos == 123241793) |
               (gene == "EML4"  & grepl("REPEAT_SEQ=T", notes) & grepl("LENGTH=7", notes) & chrom == 2  & start_pos == 42524359) |
               (gene == "EML4"  & grepl("REPEAT_SEQ=T", notes) & grepl("LENGTH=7", notes) & chrom == 2  & start_pos == 42525824) |
               (gene == "ROS1"  & grepl("REPEAT_SEQ=A", notes) & grepl("LENGTH=7", notes) & chrom == 6  & start_pos == 117658283) |
               (gene == "BRCA2" & grepl("REPEAT_SEQ=A", notes) & grepl("LENGTH=7", notes) & chrom == 13 & start_pos == 32911442) |
               (gene == "KRAS"  & grepl("REPEAT_SEQ=T", notes) & grepl("LENGTH=7", notes) & chrom == 12 & start_pos == 25368433) |
               (gene == "BRCA2" & grepl("REPEAT_SEQ=A", notes) & grepl("LENGTH=7", notes) & chrom == 13 & start_pos == 32906602)) %>%
      mutate(putative_call = "negative")

    # ghm snp output
    i_ghm <- i_ghm %>%
      filter(
        # snv
        (chrom == 7 & pos == 55259515) |
          (chrom == 7 & pos == 55249071) |
          (chrom == 7 & pos == 116412044) |
          (chrom == 17 & pos == 7579582) |
          (chrom == 7 & pos == 55242465) |

          # MSI
          (chrom == i_msi_pos$chrom[1] & pos == i_msi_pos$start_pos[1]) |
          (chrom == i_msi_neg$chrom[1] & pos == i_msi_neg$start_pos[1]) |
          (chrom == i_msi_neg$chrom[2] & pos == i_msi_neg$start_pos[2]) |
          (chrom == i_msi_neg$chrom[3] & pos == i_msi_neg$start_pos[3]) |
          (chrom == i_msi_neg$chrom[4] & pos == i_msi_neg$start_pos[4]) |
          (chrom == i_msi_neg$chrom[5] & pos == i_msi_neg$start_pos[5]) |
          (chrom == i_msi_neg$chrom[6] & pos == i_msi_neg$start_pos[6]) |
          (chrom == i_msi_neg$chrom[7] & pos == i_msi_neg$start_pos[7]) |
          (chrom == i_msi_neg$chrom[8] & pos == i_msi_neg$start_pos[8]) |
          (chrom == i_msi_neg$chrom[9] & pos == i_msi_neg$start_pos[9]) |
          (chrom == i_msi_neg$chrom[10] & pos == i_msi_neg$start_pos[10]) |

          (chrom == snv_indel_cnv_negatives$chrom[1] & pos == snv_indel_cnv_negatives$position[1]) |
          (chrom == snv_indel_cnv_negatives$chrom[2] & pos == snv_indel_cnv_negatives$position[2]) |
          (chrom == snv_indel_cnv_negatives$chrom[3] & pos == snv_indel_cnv_negatives$position[3]) |
          (chrom == snv_indel_cnv_negatives$chrom[4] & pos == snv_indel_cnv_negatives$position[4]) |
          (chrom == snv_indel_cnv_negatives$chrom[5] & pos == snv_indel_cnv_negatives$position[5]) |
          (chrom == snv_indel_cnv_negatives$chrom[6] & pos == snv_indel_cnv_negatives$position[6]) |
          (chrom == snv_indel_cnv_negatives$chrom[7] & pos == snv_indel_cnv_negatives$position[7]) |
          (chrom == snv_indel_cnv_negatives$chrom[8] & pos == snv_indel_cnv_negatives$position[8]) |
          (chrom == snv_indel_cnv_negatives$chrom[9] & pos == snv_indel_cnv_negatives$position[9]) |
          (chrom == snv_indel_cnv_negatives$chrom[10] & pos == snv_indel_cnv_negatives$position[10]) |
          (chrom == snv_indel_cnv_negatives$chrom[11] & pos == snv_indel_cnv_negatives$position[11]) |
          (chrom == snv_indel_cnv_negatives$chrom[12] & pos == snv_indel_cnv_negatives$position[12]) |
          (chrom == snv_indel_cnv_negatives$chrom[13] & pos == snv_indel_cnv_negatives$position[13]) |
          (chrom == snv_indel_cnv_negatives$chrom[14] & pos == snv_indel_cnv_negatives$position[14]) |
          (chrom == snv_indel_cnv_negatives$chrom[15] & pos == snv_indel_cnv_negatives$position[15]) |
          (chrom == snv_indel_cnv_negatives$chrom[16] & pos == snv_indel_cnv_negatives$position[16]) |
          (chrom == snv_indel_cnv_negatives$chrom[17] & pos == snv_indel_cnv_negatives$position[17]) |
          (chrom == snv_indel_cnv_negatives$chrom[18] & pos == snv_indel_cnv_negatives$position[18]) |
          (chrom == snv_indel_cnv_negatives$chrom[19] & pos == snv_indel_cnv_negatives$position[19])) %>%
      select(run_sample_id, runid, chrom, pos, non_singleton_cnt)
    
    # compile output
    i_snv_calls    <- bind_rows(i_snv_pos, i_snv_neg)
    i_indel_calls  <- bind_rows(i_indel_pos, i_indel_neg)
    i_cnv_calls    <- bind_rows(i_cnv_pos, i_cnv_neg)
    i_fusion_calls <- bind_rows(i_fusion_pos, i_fusion_neg)
    i_msi_calls <- bind_rows(i_msi_pos, i_msi_neg)
    
    output <- list(snv = i_snv_calls, 
                   indel = i_indel_calls, 
                   cnv = i_cnv_calls, 
                   fusion = i_fusion_calls,
                   msi = i_msi_calls, 
                   ghm_snp = i_ghm)
    
    return(output)
  })
  
})

# transpose the inner lists of this object so we can create the individual variant class line data sheets
# that the data analysis script expects
variant_agreement_t <- map(variant_agreement, ~transpose(.x))

flowcell_qc <- map_dfr(fc_dir, ~{
  read_tsv(file.path(.x, "autoqc_sample_qc.hdr.tsv"), col_types = col_types_autoqc) %>%
  filter(category == "flowcell") 
  })
control_qc <- map_dfr(fc_dir, ~{
  read_tsv(file.path(.x, "autoqc_sample_qc.hdr.tsv"), col_types = col_types_autoqc) %>%
  filter(category == "control")
  })
sample_qc <- map_dfr(fc_dir, ~{
  read_tsv(file.path(.x, "autoqc_sample_qc.hdr.tsv"), col_types = col_types_autoqc) %>%
  filter(category == "sample") # %>% 
  # mutate(lot = ifelse(run_sample_id %in% paste0("VCv2R", 12:22), "B", "A"))
  })

snv = c()
indel = c()
cnv = c()
fusion = c()
msi = c()
ghm_snp = c()
for (i in seq_along(variant_agreement_t)) {
  snv <- c(snv, variant_agreement_t[[i]]$snv)
  indel <- c(indel, variant_agreement_t[[i]]$indel)
  cnv <- c(cnv, variant_agreement_t[[i]]$cnv)
  fusion <- c(fusion, variant_agreement_t[[i]]$fusion)
  msi <- c(msi, variant_agreement_t[[i]]$msi)
  ghm_snp <- c(ghm_snp, variant_agreement_t[[i]]$ghm_snp)
}

snv <- do.call(rbind, snv)
indel <- do.call(rbind, indel)
cnv <- do.call(rbind, cnv)
fusion <- do.call(rbind, fusion)
msi <- do.call(plyr::rbind.fill, msi)
#msi <- do.call(rbind, c(msi, fill = TRUE))
ghm_snp <- do.call(rbind, ghm_snp)

line_data <- list('Flowcell QC' = flowcell_qc,
                  'Control QC' = control_qc, 
                  'Sample QC' = sample_qc, 
                  'SNV Calls' =   snv, 
                  'Indel Calls' =  indel, 
                  'CNV Calls' =    cnv, 
                  'Fusion Calls' = fusion, 
                  'MSI' = msi, 
                  'GHM SNP' = ghm_snp)


#line_data <- list('Flowcell QC' = flowcell_qc,
#                  'Control QC' = control_qc, 
#                  'Sample QC' = sample_qc, 
#                  'SNV Calls' =   map_dfr(variant_agreement_t, ~pluck(.x, "snv")), 
#                  'Indel Calls' =  map_dfr(variant_agreement_t, ~pluck(.x, "indel")), 
#                  'CNV Calls' =    map_dfr(variant_agreement_t, ~pluck(.x, "cnv")), 
#                  'Fusion Calls' = map_dfr(variant_agreement_t, ~pluck(.x, "fusion")), 
#                  'MSI' = map_dfr(variant_agreement_t, ~pluck(.x, "msi")), 
#                  'GHM SNP' = map_dfr(variant_agreement_t, ~pluck(.x, "ghm_snp")))


#line_data <- list('Flowcell QC' = flowcell_qc,
#                  'Control QC' = control_qc, 
#                  'Sample QC' = sample_qc, 
#                  'SNV Calls' =   map(variant_agreement_t, ~pluck(.x, "snv")) %>% bind_rows(), 
#                  'Indel Calls' =  indel <- map(variant_agreement_t, ~pluck(.x, "indel")) %>% bind_rows(), 
#                  'CNV Calls' =    cnv <- map(variant_agreement_t, ~pluck(.x, "cnv")) %>% bind_rows(), 
#                  'Fusion Calls' = fusion <- map(variant_agreement_t, ~pluck(.x, "fusion")) %>% bind_rows(),
#                  'MSI' = map(variant_agreement_t, ~pluck(.x, "msi")) %>% bind_rows(), 
#                  'GHM SNP' = map(variant_agreement_t, ~pluck(.x, "ghm_snp")) %>% bind_rows())


# write line data
openxlsx::write.xlsx(x = line_data, 
                     file = "TLFs/line_data_shanfu_wu.xlsx")
