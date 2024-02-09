# 
# step02_g360_analysis.R
# 

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("/home/swu/biostat_sns_20230207_rpt006876_g360_vcv2_performance_control_verification_report-main/")


# load libraries
library(tidyverse)
library(openxlsx)
library(beanplot)

source("Code/myblend.R")

# load study data
snv    <- openxlsx::read.xlsx("TLFs/line_data_shanfu_wu.xlsx", sheet = "SNV Calls") %>%
  mutate(id1 = gsub("^(.+)_(.+)_(.+)_(.+)", "\\1", run_sample_id),
         tp = gsub("^(.+)_(.+)_(.+)_(.+)", "\\2", run_sample_id))
indel  <- openxlsx::read.xlsx("TLFs/line_data_shanfu_wu.xlsx", sheet = "Indel Calls") %>%
  mutate(id1 = gsub("^(.+)_(.+)_(.+)_(.+)", "\\1", run_sample_id),
         tp = gsub("^(.+)_(.+)_(.+)_(.+)", "\\2", run_sample_id))
cnv    <- openxlsx::read.xlsx("TLFs/line_data_shanfu_wu.xlsx", sheet = "CNV Calls") %>%
  mutate(id1 = gsub("^(.+)_(.+)_(.+)_(.+)", "\\1", run_sample_id),
         tp = gsub("^(.+)_(.+)_(.+)_(.+)", "\\2", run_sample_id))
fusion <- openxlsx::read.xlsx("TLFs/line_data_shanfu_wu.xlsx", sheet = "Fusion Calls") %>%
  mutate(id1 = gsub("^(.+)_(.+)_(.+)_(.+)", "\\1", run_sample_id),
         tp = gsub("^(.+)_(.+)_(.+)_(.+)", "\\2", run_sample_id))
msi    <- openxlsx::read.xlsx("TLFs/line_data_shanfu_wu.xlsx", sheet = "MSI") %>%
  mutate(id1 = gsub("^(.+)_(.+)_(.+)_(.+)", "\\1", run_sample_id),
         tp = gsub("^(.+)_(.+)_(.+)_(.+)", "\\2", run_sample_id))
ghm    <- openxlsx::read.xlsx("TLFs/line_data_shanfu_wu.xlsx", sheet = "GHM SNP") %>%
  mutate(id1 = gsub("^(.+)_(.+)_(.+)_(.+)", "\\1", run_sample_id),
         tp = gsub("^(.+)_(.+)_(.+)_(.+)", "\\2", run_sample_id))
snv_indel_msi <- snv %>% 
  distinct(chrom, position, gene, mut_nt, putative_call) %>%
  mutate(variant_class = "snv") %>%
  bind_rows(indel %>% 
              distinct(chrom, position, gene, mut_nt, putative_call) %>%
              mutate(variant_class = "indel")) %>%
  bind_rows(msi %>% 
              distinct(chrom, position = start_pos, gene, putative_call) %>%
              mutate(variant_class = "msi"))

ghm <- ghm %>%
  left_join(snv_indel_msi, by = c("chrom", "pos" = "position"))


# ---------------------------------------------------------
# -----------------------------------------------------
# STEP 1: Primary analysis (PPA and NPA)
# -----------------------------------------------------
# ---------------------------------------------------------

# ------------------------------------------
# SNV Analysis
# ------------------------------------------
snv_summary_pos <- snv %>% 
  filter(putative_call == "positive") %>%
  group_by(id1, tp) %>%
  summarise(x = sum(call == 1), 
            n = n(), 
            ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)")) %>%
  mutate(var_class = "snv")

snv_summary_neg <- snv %>% 
  filter(putative_call == "negative") %>%
  group_by(id1, tp) %>%
  summarise(x = sum(call == 0 | is.na(call)), 
            n = n(), 
            ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)")) %>%
  mutate(var_class = "snv")

# ------------------------------------------
# Indel Analysis
# ------------------------------------------

indel_summary_pos <- indel %>% 
  filter(putative_call == "positive") %>%
  group_by(id1, tp) %>%
  summarise(x = sum(call == 1), 
            n = n(), 
            ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)")) %>%
  mutate(var_class = "indel")

indel_summary_neg <- indel %>% 
  filter(putative_call == "negative") %>%
  group_by(id1, tp) %>%
  summarise(x = sum(call == 0 | is.na(call)), 
            n = n(), 
            ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)")) %>%
  mutate(var_class = "indel")

# ------------------------------------------
# CNV Analysis
# ------------------------------------------

cnv_summary_pos <- cnv %>% 
  filter(putative_call == "positive") %>%
  group_by(id1, tp) %>%
  summarise(x = sum(call == 2), 
            n = n(), 
            ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)")) %>%
  mutate(var_class = "cnv")

cnv_summary_neg <- cnv %>% 
  filter(putative_call == "negative") %>%
  group_by(id1, tp) %>%
  summarise(x = sum(call == 0 | is.na(call)), 
            n = n(), 
            ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)")) %>%
  mutate(var_class = "cnv")

# ------------------------------------------
# Fusion Analysis
# ------------------------------------------

fusion_summary_pos <- fusion %>% 
  filter(putative_call == "positive") %>%
  group_by(id1, tp) %>%
  summarise(x = sum(call == 1), 
            n = n(), 
            ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)")) %>%
  mutate(var_class = "fusion")

fusion_summary_neg <- fusion %>% 
  filter(putative_call == "negative") %>%
  group_by(id1, tp) %>%
  summarise(x = sum(call == 0 | is.na(call)), 
            n = n(), 
            ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
                         sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)")) %>%
  mutate(var_class = "fusion")

# ------------------------------------------
# Acceptance Criteria Assessment
# ------------------------------------------
study_ppa <- bind_rows(snv_summary_pos, indel_summary_pos, cnv_summary_pos, fusion_summary_pos) %>%
  group_by(id1, tp) %>%
  summarise(x = sum(x),
            n = sum(n),
            PointEst = Hmisc::binconf(x = sum(x), n = sum(n), method = "exact")[,1],
            Lower = Hmisc::binconf(x = sum(x), n = sum(n), method = "exact")[,2],
            Upper = Hmisc::binconf(x = sum(x), n = sum(n), method = "exact")[,3]) %>%
  ungroup() %>%
  mutate(metric = "PPA") %>%
  select(metric, id1, tp, x, n, PointEst, Lower, Upper)

study_npa <- bind_rows(snv_summary_neg, indel_summary_neg, cnv_summary_neg, fusion_summary_neg) %>%
  group_by(id1, tp) %>%
  summarise(x = sum(x),
            n = sum(n),
            PointEst = Hmisc::binconf(x = sum(x), n = sum(n), method = "exact")[,1],
            Lower = Hmisc::binconf(x = sum(x), n = sum(n), method = "exact")[,2],
            Upper = Hmisc::binconf(x = sum(x), n = sum(n), method = "exact")[,3]) %>%
  ungroup() %>%
  mutate(metric = "NPA") %>%
  select(metric, id1, tp, x, n, PointEst, Lower, Upper)

# write summary table for primary analysis
write_csv(x = rbind.data.frame(study_ppa, study_npa) %>%
            mutate(PointEst = round(100*PointEst, 1), 
                   Lower = round(100*Lower, 2), 
                   Upper = round(100*Upper, 2)),
          "TLFs/NaOH_Test_D30_agreement_summary.csv")


# # ---------------------------------------------------------
# # -----------------------------------------------------
# # STEP 2: Exploratory Analysis
# # -----------------------------------------------------
# # ---------------------------------------------------------
# 
# # ------------------------
# # Step 2a: MSI Agreement
# # ------------------------
# 
# msi_summary_pos <- msi %>% 
#   filter(putative_call == "positive") %>%
#   summarise(x = sum(delta_aic > 5), 
#             n = n(), 
#             agreement = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)"))
# 
# msi_summary_neg <- msi %>% 
#   filter(putative_call == "negative") %>%
#   summarise(x = sum(delta_aic <= 5), 
#             n = n(), 
#             agreement = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)"))
# 
# msi_summaries <- bind_rows(msi_summary_pos, 
#                            msi_summary_neg) %>%
#   mutate(putative_call = c("positive", "negative")) %>%
#   select(putative_call, x, n, agreement)
# 
# readr::write_csv(msi_summaries, path = "TLFs/exploratory_msi_site_agreement.csv")
# 
# if(F){
#   ##################################################
#   # Not needed since positive MSI site was updated
#   ##################################################
# 
#   # -----------------------
#   # MSI Positive deltaAIC
#   # -----------------------
#   
#   png("TLFs/beanplot_msi_putative_positive_deltaAIC.png",
#       width = 5, height = 5, units = "in", pointsize = 10, res = 600)
#   par(mar = c(1,4,4,1))
#   par(mgp = c(2,0.6,0))
#   beanplot(msi$delta_aic[msi$putative_call == "positive"], 
#            col = list(c(myblend(c(1,1), "hotpink", "gray90"), "black", "black", "black"), 
#                       c("yellow", "black", "black", "black")), 
#            ll = 0.025, what = c(0,1,0,1), 
#            axes = FALSE)
#   axis(2, las = 2)
#   title(ylab = "deltaAIC")
#   title(main = "Distribution of deltaAIC\nfor MSI Putative Positive Sites")
#   abline(h = 5, col = myblend(c(1,1), "turquoise4", "gray90"), lty = "dashed")
#   text(x = 0.5, y = 5.5, labels = "MSI Positive", pos = 4, srt = 90)
#   text(x = 0.5, y = 4.5, labels = "MSI Negative", pos = 4, srt = 270)
#   dev.off()
# } 
# 
# # ---------------------------
# # Step 2b: Coverage Summary
# # ---------------------------
# 
# # summarise data from ghm_snp
# coverage_table_by_variant <- ghm %>%
#   mutate(variant_class = factor(variant_class, 
#                                 levels = c("snv", "indel", "msi"), 
#                                 labels = c("SNV", "Indel", "MSI"))) %>%
#   distinct(run_sample_id, variant_class, putative_call, gene, chrom, pos, .keep_all = TRUE) %>%
#   group_by(variant_class, putative_call, gene, chrom, position = pos) %>%
#   summarise(mean_coverage = round(mean(non_singleton_cnt), 1), 
#             median_coverage = round(median(non_singleton_cnt), 1), 
#             sd_coverage = round(sd(non_singleton_cnt), 1))
# 
# # write coverage summary table for snv, indel, msi
# write_csv(coverage_table_by_variant, 
#           path = "TLFs/exploratory_coverage_summaries_for_snv_indel_msi_sites.csv")
# 
# # print summary table for coverages > 1000
# ghm %>% 
#   distinct(run_sample_id, variant_class, putative_call, gene, chrom, pos, .keep_all = TRUE) %>%
#   summarise(prop_greater_1000 = mean(non_singleton_cnt > 1000), 
#             n_greater_1000 = sum(non_singleton_cnt > 1000), 
#             total = n())
# 
# # write coverage summary table for fusions
# fusion %>% 
#   group_by(gene_a, gene_b) %>%
#   summarise(mean_cov_a = mean(molecule_coverage_a), 
#             mean_cov_b = mean(molecule_coverage_b), 
#             median_cov_a = median(molecule_coverage_a), 
#             median_cov_b = median(molecule_coverage_b), 
#             sd_cov_a = sd(molecule_coverage_a), 
#             sd_cov_b = sd(molecule_coverage_b)) %>%
#   write_csv(path = "TLFs/exploratory_coverage_summaries_for_fusions.csv")
# 
# 
# # ---------------------------
# # Step 2c: Allele frequency
# # ---------------------------
# 
# # write allele frequency summary table for positive SNV, indel, cnv, and fusion sites
# snv %>%
#   mutate(class = "SNV") %>%
#   filter(putative_call == "positive") %>%
#   group_by(class, gene, chrom, position) %>%
#   summarise(mean_maf = round(mean(percentage), 1), 
#             median_maf = round(median(percentage), 1), 
#             sd_maf = round(sd(percentage), 1)) %>%
#   bind_rows(indel %>%
#               mutate(class = "Indel") %>%
#               filter(putative_call == "positive") %>%
#               group_by(class, gene, chrom, position) %>%
#               summarise(mean_maf = round(mean(percentage), 1), 
#                         median_maf = round(median(percentage), 1), 
#                         sd_maf = round(sd(percentage), 1))) %>%
#   bind_rows(cnv %>%
#               mutate(class = "CNV") %>%
#               filter(putative_call == "positive") %>%
#               group_by(class, gene) %>%
#               summarise(mean_cn = round(mean(copy_number), 2), 
#                         median_cn = round(median(copy_number), 2), 
#                         sd_cn = round(sd(copy_number), 2))) %>%
#   bind_rows(fusion %>%
#               mutate(class = "Fusion") %>%
#               filter(putative_call == "positive") %>%
#               group_by(class, gene_a, gene_b) %>%
#               summarise(mean_maf = round(mean(percentage), 1), 
#                         median_maf = round(median(percentage), 1), 
#                         sd_maf = round(sd(percentage), 1))) %>%
#   write_csv(path = "TLFs/exploratory_allele_frequecy_summaries_for_snv_indel_fusion_cnv.csv")
# 
# 
# # --------------------------------
# # Step 2c: Repeat Steps 2a-c
# #  stratified by GHM material lot
# # --------------------------------
# 
# # ------------------------
# # Agreement rates by lot
# # ------------------------
# 
# snv_summary_pos_lot <- snv %>% 
#   filter(putative_call == "positive") %>%
#   group_by(lot) %>%
#   summarise(x = sum(call == 1), 
#             n = n(), 
#             ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)"))
# 
# snv_summary_neg_lot <- snv %>% 
#   filter(putative_call == "negative") %>%
#   group_by(lot) %>%
#   summarise(x = sum(call == 0 | is.na(call)), 
#             n = n(), 
#             ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)"))
# 
# # ------------------------------------------
# # Indel Analysis
# # ------------------------------------------
# 
# indel_summary_pos_lot <- indel %>% 
#   filter(putative_call == "positive") %>%
#   group_by(lot) %>%
#   summarise(x = sum(call == 1), 
#             n = n(), 
#             ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)"))
# 
# indel_summary_neg_lot <- indel %>% 
#   filter(putative_call == "negative") %>%
#   group_by(lot) %>%
#   summarise(x = sum(call == 0 | is.na(call)), 
#             n = n(), 
#             ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)"))
# 
# # ------------------------------------------
# # CNV Analysis
# # ------------------------------------------
# 
# cnv_summary_pos_lot <- cnv %>% 
#   filter(putative_call == "positive") %>%
#   group_by(lot) %>%
#   summarise(x = sum(call == 2), 
#             n = n(), 
#             ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)"))
# 
# cnv_summary_neg_lot <- cnv %>% 
#   filter(putative_call == "negative") %>%
#   group_by(lot) %>%
#   summarise(x = sum(call == 0 | is.na(call)), 
#             n = n(), 
#             ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)"))
# 
# # ------------------------------------------
# # Fusion Analysis
# # ------------------------------------------
# 
# fusion_summary_pos_lot <- fusion %>% 
#   filter(putative_call == "positive") %>%
#   group_by(lot) %>%
#   summarise(x = sum(call == 1), 
#             n = n(), 
#             ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)"))
# 
# fusion_summary_neg_lot <- fusion %>% 
#   filter(putative_call == "negative") %>%
#   group_by(lot) %>%
#   summarise(x = sum(call == 0 | is.na(call)), 
#             n = n(), 
#             ppa = paste0(sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"PointEst"]), "% (",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Lower"]), "%, ",
#                          sprintf("%.1f", 100*Hmisc::binconf(x, n, method = "exact")[1,"Upper"]), "%)"))
# 
# # ------------------------------------------
# # Compile output for exploratory agreements
# # ------------------------------------------
# 
# study_ppa_A <- Hmisc::binconf(x = sum(snv_summary_pos_lot$x[snv_summary_pos_lot$lot == "A"], 
#                                     indel_summary_pos_lot$x[indel_summary_pos_lot$lot == "A"], 
#                                     cnv_summary_pos_lot$x[cnv_summary_pos_lot$lot == "A"],
#                                     fusion_summary_pos_lot$x[fusion_summary_pos_lot$lot == "A"]), 
#                             n = sum(snv_summary_pos_lot$n[snv_summary_pos_lot$lot == "A"], 
#                                     indel_summary_pos_lot$n[indel_summary_pos_lot$lot == "A"], 
#                                     cnv_summary_pos_lot$n[cnv_summary_pos_lot$lot == "A"],
#                                     fusion_summary_pos_lot$n[fusion_summary_pos_lot$lot == "A"]), 
#                             method = "exact") %>% 
#   as.data.frame() %>%
#   mutate(x = sum(snv_summary_pos_lot$x[snv_summary_pos_lot$lot == "A"], 
#                   indel_summary_pos_lot$x[indel_summary_pos_lot$lot == "A"], 
#                   cnv_summary_pos_lot$x[cnv_summary_pos_lot$lot == "A"],
#                   fusion_summary_pos_lot$x[fusion_summary_pos_lot$lot == "A"]), 
#           n = sum(snv_summary_pos_lot$n[snv_summary_pos_lot$lot == "A"], 
#                   indel_summary_pos_lot$n[indel_summary_pos_lot$lot == "A"], 
#                   cnv_summary_pos_lot$n[cnv_summary_pos_lot$lot == "A"],
#                   fusion_summary_pos_lot$n[fusion_summary_pos_lot$lot == "A"]), 
#          Lot = "A") %>%
#   mutate(metric = "PPA") %>% 
#   select(metric, x, n, PointEst, Lower, Upper, Lot)
# 
# study_ppa_B <- Hmisc::binconf(x = sum(snv_summary_pos_lot$x[snv_summary_pos_lot$lot == "B"], 
#                                       indel_summary_pos_lot$x[indel_summary_pos_lot$lot == "B"], 
#                                       cnv_summary_pos_lot$x[cnv_summary_pos_lot$lot == "B"],
#                                       fusion_summary_pos_lot$x[fusion_summary_pos_lot$lot == "B"]), 
#                               n = sum(snv_summary_pos_lot$n[snv_summary_pos_lot$lot == "B"], 
#                                       indel_summary_pos_lot$n[indel_summary_pos_lot$lot == "B"], 
#                                       cnv_summary_pos_lot$n[cnv_summary_pos_lot$lot == "B"],
#                                       fusion_summary_pos_lot$n[fusion_summary_pos_lot$lot == "B"]), 
#                               method = "exact") %>% 
#   as.data.frame() %>%
#   mutate(x = sum(snv_summary_pos_lot$x[snv_summary_pos_lot$lot == "B"], 
#                  indel_summary_pos_lot$x[indel_summary_pos_lot$lot == "B"], 
#                  cnv_summary_pos_lot$x[cnv_summary_pos_lot$lot == "B"],
#                  fusion_summary_pos_lot$x[fusion_summary_pos_lot$lot == "B"]), 
#          n = sum(snv_summary_pos_lot$n[snv_summary_pos_lot$lot == "B"], 
#                  indel_summary_pos_lot$n[indel_summary_pos_lot$lot == "B"], 
#                  cnv_summary_pos_lot$n[cnv_summary_pos_lot$lot == "B"],
#                  fusion_summary_pos_lot$n[fusion_summary_pos_lot$lot == "B"]), 
#          Lot = "B") %>%
#   mutate(metric = "PPA") %>% 
#   select(metric, x, n, PointEst, Lower, Upper, Lot)
# 
# study_npa_A <- Hmisc::binconf(x = sum(snv_summary_neg_lot$x[snv_summary_neg_lot$lot == "A"], 
#                                       indel_summary_neg_lot$x[indel_summary_neg_lot$lot == "A"], 
#                                       cnv_summary_neg_lot$x[cnv_summary_neg_lot$lot == "A"],
#                                       fusion_summary_neg_lot$x[fusion_summary_neg_lot$lot == "A"]), 
#                               n = sum(snv_summary_neg_lot$n[snv_summary_neg_lot$lot == "A"], 
#                                       indel_summary_neg_lot$n[indel_summary_neg_lot$lot == "A"], 
#                                       cnv_summary_neg_lot$n[cnv_summary_neg_lot$lot == "A"],
#                                       fusion_summary_neg_lot$n[fusion_summary_neg_lot$lot == "A"]), 
#                               method = "exact") %>% 
#   as.data.frame() %>%
#   mutate(x = sum(snv_summary_neg_lot$x[snv_summary_neg_lot$lot == "A"], 
#                  indel_summary_neg_lot$x[indel_summary_neg_lot$lot == "A"], 
#                  cnv_summary_neg_lot$x[cnv_summary_neg_lot$lot == "A"],
#                  fusion_summary_neg_lot$x[fusion_summary_neg_lot$lot == "A"]), 
#          n = sum(snv_summary_neg_lot$n[snv_summary_neg_lot$lot == "A"], 
#                  indel_summary_neg_lot$n[indel_summary_neg_lot$lot == "A"], 
#                  cnv_summary_neg_lot$n[cnv_summary_neg_lot$lot == "A"],
#                  fusion_summary_neg_lot$n[fusion_summary_neg_lot$lot == "A"]), 
#          Lot = "A") %>%
#   mutate(metric = "NPA") %>% 
#   select(metric, x, n, PointEst, Lower, Upper, Lot)
# 
# study_npa_B <- Hmisc::binconf(x = sum(snv_summary_neg_lot$x[snv_summary_neg_lot$lot == "B"], 
#                                       indel_summary_neg_lot$x[indel_summary_neg_lot$lot == "B"], 
#                                       cnv_summary_neg_lot$x[cnv_summary_neg_lot$lot == "B"],
#                                       fusion_summary_neg_lot$x[fusion_summary_neg_lot$lot == "B"]), 
#                               n = sum(snv_summary_neg_lot$n[snv_summary_neg_lot$lot == "B"], 
#                                       indel_summary_neg_lot$n[indel_summary_neg_lot$lot == "B"], 
#                                       cnv_summary_neg_lot$n[cnv_summary_neg_lot$lot == "B"],
#                                       fusion_summary_neg_lot$n[fusion_summary_neg_lot$lot == "B"]), 
#                               method = "exact") %>% 
#   as.data.frame() %>%
#   mutate(x = sum(snv_summary_neg_lot$x[snv_summary_neg_lot$lot == "B"], 
#                  indel_summary_neg_lot$x[indel_summary_neg_lot$lot == "B"], 
#                  cnv_summary_neg_lot$x[cnv_summary_neg_lot$lot == "B"],
#                  fusion_summary_neg_lot$x[fusion_summary_neg_lot$lot == "B"]), 
#          n = sum(snv_summary_neg_lot$n[snv_summary_neg_lot$lot == "B"], 
#                  indel_summary_neg_lot$n[indel_summary_neg_lot$lot == "B"], 
#                  cnv_summary_neg_lot$n[cnv_summary_neg_lot$lot == "B"],
#                  fusion_summary_neg_lot$n[fusion_summary_neg_lot$lot == "B"]), 
#          Lot = "B") %>%
#   mutate(metric = "NPA") %>% 
#   select(metric, x, n, PointEst, Lower, Upper, Lot)
# 
# write_csv(x = rbind.data.frame(study_ppa_A, study_ppa_B, study_npa_B, study_npa_B) %>%
#             select(Lot, metric, x, n, PointEst, Lower, Upper) %>%
#             mutate(PointEst = paste0(round(100*PointEst, 1), "%"), 
#                    Lower = paste0(round(100*Lower, 2), "%"), 
#                    Upper = paste0(round(100*Upper, 2), "%")),
#           path = "TLFs/g360_agreement_summary_by_lot.csv")
# 
# 
# # ------------------------
# # Coverage by lot
# # ------------------------
# 
# coverage_table_by_variant_lot <- ghm %>%
#   mutate(variant_class = factor(variant_class, 
#                                 levels = c("snv", "indel", "msi"), 
#                                 labels = c("SNV", "Indel", "MSI"))) %>%
#   distinct(run_sample_id, variant_class, putative_call, lot, gene, chrom, pos, .keep_all = TRUE) %>%
#   group_by(variant_class, lot, putative_call, gene, chrom, position = pos) %>%
#   summarise(mean_coverage = round(mean(non_singleton_cnt), 1), 
#             median_coverage = round(median(non_singleton_cnt), 1), 
#             sd_coverage = round(sd(non_singleton_cnt), 1)) %>%
#   arrange(variant_class, putative_call, gene, chrom, position, lot)
# 
# write_csv(coverage_table_by_variant_lot, 
#           path = "TLFs/exploratory_coverage_summaries_for_snv_indel_msi_sites_by_lot.csv")
# 
# ghm %>% 
#   distinct(run_sample_id, variant_class, putative_call, gene, chrom, pos, lot, .keep_all = TRUE) %>%
#   group_by(lot) %>%
#   summarise(prop_greater_1000 = mean(non_singleton_cnt > 1000), 
#             n_greater_1000 = sum(non_singleton_cnt > 1000), 
#             total = n())
# 
# fusion %>% 
#   group_by(gene_a, gene_b, lot) %>%
#   summarise(mean_cov_a = mean(molecule_coverage_a), 
#             mean_cov_b = mean(molecule_coverage_b), 
#             median_cov_a = median(molecule_coverage_a), 
#             median_cov_b = median(molecule_coverage_b), 
#             sd_cov_a = sd(molecule_coverage_a), 
#             sd_cov_b = sd(molecule_coverage_b)) %>%
#   write_csv(path = "TLFs/exploratory_coverage_summaries_for_fusions_by_lot.csv")
# 
# 
# # ------------------------
# # Allele freq by lot
# # ------------------------
# 
# snv %>%
#   mutate(class = "SNV") %>%
#   filter(putative_call == "positive") %>%
#   group_by(class, gene, chrom, position, lot) %>%
#   summarise(mean_maf = round(mean(percentage), 1), 
#             median_maf = round(median(percentage), 1), 
#             sd_maf = round(sd(percentage), 1)) %>%
#   bind_rows(indel %>%
#               mutate(class = "Indel") %>%
#               filter(putative_call == "positive") %>%
#               group_by(class, gene, chrom, position, lot) %>%
#               summarise(mean_maf = round(mean(percentage), 1), 
#                         median_maf = round(median(percentage), 1), 
#                         sd_maf = round(sd(percentage), 1))) %>%
#   bind_rows(cnv %>%
#               mutate(class = "CNV") %>%
#               filter(putative_call == "positive") %>%
#               group_by(class, gene, lot) %>%
#               summarise(mean_cn = round(mean(copy_number), 2), 
#                         median_cn = round(median(copy_number), 2), 
#                         sd_cn = round(sd(copy_number), 2))) %>%
#   bind_rows(fusion %>%
#               mutate(class = "Fusion") %>%
#               filter(putative_call == "positive") %>%
#               group_by(class, gene_a, gene_b, lot) %>%
#               summarise(mean_maf = round(mean(percentage), 1), 
#                         median_maf = round(median(percentage), 1), 
#                         sd_maf = round(sd(percentage), 1))) %>%
#   write_csv(path = "TLFs/exploratory_allele_frequecy_summaries_for_snv_indel_fusion_cnv_by_lot.csv")
# 
