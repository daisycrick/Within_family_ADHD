
rm(list = ls())
setwd ("~/OneDrive - University College London/Documents/01_projects/01_ADHD_brisbane_paper") # Set the file path to your local directory

#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")
library(dplyr)
library(TwoSampleMR)
library(purrr)
library(ieugwasr)


#---- Load and Prepare Exposure Data ----#
ADHD <- read_exposure_data(
  filename = "27_snps_condensed.txt",
  sep = " ",
  snp_col = "SNP", 
  beta_col = "logOR", 
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2", 
  eaf_col = "EAF",
  pval_col = "PVALUE"
)

#---- Function to run MR for one outcome ----#
run_mr_pipeline <- function(outcome_id, exposure_data, outcome_sample_size, exposure_sample_size = 487410, output_prefix = "./population/reduced_population") {
  # Extract and harmonise
  outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = outcome_id, proxies = TRUE)
  


 snps_to_remove_by_outcome <- list(
   "ukb-b-5779" = c("rs114142727", "rs73145587", "rs76284431"),
  "ieu-b-109" = c("rs114142727", "rs73145587", "rs76284431"),
   "ieu-b-110" = c("rs114142727", "rs73145587", "rs76284431"),
   "ieu-b-111" = c("rs114142727", "rs73145587", "rs76284431"),
  "ieu-b-142" = c("rs114142727", "rs73145587", "rs76284431"),
   "ieu-a-1001" = c("rs114142727", "rs73145587", "rs76284431"),
   "ukb-b-17422" = c("rs114142727", "rs73145587", "rs76284431"),
  "ebi-a-GCST90029028" = c("rs114142727", "rs76284431", "rs73145587"),
  "ieu-b-40" = c("rs76284431"),
  "ieu-a-1000" = c("rs114142727","rs73145587","rs76284431"),
  "ieu-a-74" = c("rs114142727","rs73145587","rs76284431"))
  
  if (outcome_id %in% names(snps_to_remove_by_outcome)) {
    outcome_data <- outcome_data[
      !outcome_data$SNP %in% snps_to_remove_by_outcome[[outcome_id]],
    ]
  }
  
 dat <- harmonise_data(exposure_data, outcome_data, action = 2)
  
#  # Filter palindromic SNPs
  keep_snps <- subset(dat, mr_keep == TRUE)
  ADHD_filtered <- merge(exposure_data, keep_snps[, c("SNP", "mr_keep")], by = "SNP")
 ADHD_filtered <- subset(ADHD_filtered, mr_keep == TRUE)
 
 snps_before <- ADHD_filtered$SNP
 
#  # Clump
  ADHD_clumped <- clump_data(ADHD_filtered, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR")
 dat <- harmonise_data(ADHD_clumped, outcome_data, action = 2)
 
 # Identify removed SNPs
 removed_snps <- setdiff(snps_before, ADHD_clumped$SNP)

 print(removed_snps)
#  
#  # Add sample sizes
  dat$samplesize.exposure <- exposure_sample_size
  dat$samplesize.outcome <- outcome_sample_size
 # write.csv(dat, file = paste0(output_prefix, "/harmonized/harmonised_", outcome_id, ".csv"), row.names = FALSE)
#  # Run MR
  mr_result <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  mr_result$lower <- mr_result$b - 1.96 * mr_result$se
 mr_result$upper <- mr_result$b + 1.96 * mr_result$se
 
 mr_result
#  write.csv(mr_result, file = paste0(output_prefix, "/mr_", outcome_id, "_presteiger.csv"), row.names = FALSE)
  
  # Sensitivity
 #write.csv(mr_heterogeneity(dat), paste0(output_prefix, "/pre_ste_het_", outcome_id, ".csv"), row.names = FALSE)
 #write.csv(mr_pleiotropy_test(dat), paste0(output_prefix, "/pre_ste_pleio_sibling_", outcome_id, ".csv"), row.names = FALSE)
 set.seed(12345)
 presso <- run_mr_presso(dat, NbDistribution = 10000, SignifThreshold = 0.05)
 print(presso)
 
 print(dat$SNP[c(1, 4, 7, 11, 16)])
 print(dat$SNP[c(3, 12, 14, 18, 20)])
 
  
  # Print the names of the SNPs driving the heterogeneity
  
  # Add R^2
#  dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
#  dat$r.outcome <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)
  
  # Steiger filtering
#  dat$MAF <- 1 - dat$eaf.exposure
#  dat$r2 <- 2 * dat$MAF * (1 - dat$MAF) * (dat$beta.exposure^2)
#  dat$F.stat <- (dat$r2 * (dat$samplesize.exposure - 2)) / (1 - dat$r2)
#  dat <- steiger_filtering(dat)
#  write.csv(dat, file = paste0(output_prefix, "poststeiger/reduced/populations/harmonised_", outcome_id, ".csv"), row.names = FALSE)
#  dat <- subset(dat, steiger_dir == TRUE)
#  write.csv(dat, file = paste0(output_prefix, "poststeiger/reduced/populations/mr_", outcome_id, "_presteiger.csv"), row.names = FALSE)
  
  
  
  # Run MR
#  mr_result <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
#  mr_result$lower <- mr_result$b - 1.96 * mr_result$se
#  mr_result$upper <- mr_result$b + 1.96 * mr_result$se
# write.csv(mr_result, file = paste0(output_prefix, "poststeiger/reduced/populations/mr_results/mr_", outcome_id, "_poststeiger.csv"), row.names = FALSE)
  
  # Sensitivity
#  write.csv(mr_heterogeneity(dat), paste0(output_prefix, "poststeiger/reduced/populations/heterogeneity/het_", outcome_id, ".csv"), row.names = FALSE)
#  write.csv(mr_pleiotropy_test(dat), paste0(output_prefix, "poststeiger/reduced/populations/pleiotropy/pleio_", outcome_id, ".csv"), row.names = FALSE)
  #presso <- run_mr_presso(dat, NbDistribution = 10000, SignifThreshold = 0.05)
  #print(presso)
}

#---- Run the pipeline for one or multiple outcomes ----#

# OR if you want to run across multiple outcomes:

outcomes <- c("ieu-b-142", "ieu-a-1001", "ukb-b-17422", "ebi-a-GCST90029028", 
              "ieu-b-40", "ieu-a-1000", "ieu-a-74",
              "ieu-b-40", "ieu-a-74", "ieu-b-109", "ieu-b-110",
              "ieu-b-111", "ebi-a-GCST90029070", "ebi-a-GCST90014006", "ukb-b-5779",
              "ieu-b-142", "ieu-a-1001", "ebi-a-GCST90000050", "ieu-b-5106", 
              "ukb-b-17422", "ebi-a-GCST90029028", "ieu-a-1000", "ebi-a-GCST006250")

samplesizes <- c(249752, 293723, 143819, 374323,
                 681375, 161460, 124591,
                 681375, 124591, 403943, 440546,
                 441016, 575531, 389889, 462246,
                 249752, 293723, 542901, 252514,
                 143819, 374323, 161460, 269867)


for (i in seq_along(outcomes)) {
  run_mr_pipeline(outcomes[i], ADHD, outcome_sample_size = samplesizes[i])
}

