rm(list = ls())
setwd ("/Users/daisy/OneDrive - University College London/Documents/01_projects/01_ADHD_brisbane_paper/") # Set the file path to your local directory

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
run_mr_pipeline <- function(outcome_id, exposure_data, outcome_sample_size, exposure_sample_size = 487410, output_prefix = "./siblings/reduced_sample") {
  # Extract and harmonise
  outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = outcome_id, proxies = TRUE)
 
  
  #snps_to_remove_by_outcome <- list(
  #  "ieu-b-4815" = c("rs1162202")
 #)
  

  #if (outcome_id %in% names(snps_to_remove_by_outcome)) {
  #  outcome_data <- outcome_data[
  #    !outcome_data$SNP %in% snps_to_remove_by_outcome[[outcome_id]],
  #  ]
 # }
  
   dat <- harmonise_data(exposure_data, outcome_data, action = 2)
  print(dat)
  # Filter palindromic SNPs
   keep_snps <- subset(dat, mr_keep == TRUE)
   ADHD_filtered <- merge(exposure_data, keep_snps[, c("SNP", "mr_keep")], by = "SNP")
    ADHD_filtered <- subset(ADHD_filtered, mr_keep == TRUE)
snps_before <- ADHD_filtered$SNP

  # Clump
  ADHD_clumped <- clump_data(ADHD_filtered, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR")
  dat <- harmonise_data(ADHD_clumped, outcome_data, action = 2)
  
  # Identify removed SNPs
  removed_snps <- setdiff(snps_before, ADHD_clumped$SNP)
  
  print(removed_snps)
  #  
  
  # Add sample sizes
  dat$samplesize.exposure <- exposure_sample_size
  dat$samplesize.outcome <- outcome_sample_size
  write.csv(dat, file = paste0(output_prefix, "/harmonized/harmonised_", outcome_id, ".csv"), row.names = FALSE)
  mr_result <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
  # Run MR
  mr_result$lower <- mr_result$b - 1.96 * mr_result$se
  mr_result$upper <- mr_result$b + 1.96 * mr_result$se
 write.csv(mr_result, file = paste0(output_prefix, "/mr_results/mr_sibling_", outcome_id, "_presteiger.csv"), row.names = FALSE)
  
  # Sensitivity
  #write.csv(mr_pleiotropy_test(dat), paste0(output_prefix, "/pre_ste_pleio_sibling_", outcome_id, ".csv"), row.names = FALSE)
#write.csv(mr_heterogeneity(dat), paste0(output_prefix, "/het_pleio/pre_ste_het_sibling_", outcome_id, ".csv"), row.names = FALSE)
  presso <- run_mr_presso(dat, NbDistribution = 10000, SignifThreshold = 0.05)
  print(presso)
  
  
  # Steiger filtering
# dat$MAF <- 1 - dat$eaf.exposure
#  dat$r2 <- 2 * dat$MAF * (1 - dat$MAF) * (dat$beta.exposure^2)
#  dat$F.stat <- (dat$r2 * (dat$samplesize.exposure - 2)) / (1 - dat$r2)
 # dat <- steiger_filtering(dat)
#  write.csv(dat, file = paste0(output_prefix, "poststeiger/reduced/siblings/harmonised_", outcome_id, ".csv"), row.names = FALSE)
#  dat <- subset(dat, steiger_dir == TRUE)
  
  # Run MR
#  mr_result <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
 # mr_result$lower <- mr_result$b - 1.96 * mr_result$se
#  mr_result$upper <- mr_result$b + 1.96 * mr_result$se
##  write.csv(mr_result, file = paste0(output_prefix, "poststeiger/reduced/siblings/mr_results/mr_sibling_", outcome_id, "_poststeiger.csv"), row.names = FALSE)
  
  # Sensitivity
 # write.csv(mr_heterogeneity(dat), paste0(output_prefix, "poststeiger/reduced/siblings/heterogeneity/het_sibling_", outcome_id, ".csv"), row.names = FALSE)
#  write.csv(mr_pleiotropy_test(dat), paste0(output_prefix, "poststeiger/reduced/siblings/pleiotropy/pleio_sibling_", outcome_id, ".csv"), row.names = FALSE)
 #  presso <- run_mr_presso(dat, NbDistribution = 10000, SignifThreshold = 0.05)
#    print(presso)
}

#---- Run the pipeline for  multiple outcomes ----#


outcomes <- c("ieu-b-4819",
              "ieu-b-4829", "ieu-b-4843", "ieu-b-4845", "ieu-b-4815"
              "ieu-b-4849", "ieu-b-4763", "ieu-b-4841", "ieu-b-4833", 
              "ieu-b-4825", "ieu-b-4835", , "ieu-b-4821",
              "ieu-b-4823", "ieu-b-4848", "ieu-b-4839", "ieu-b-4837")
              
               # Example outcome IDs



samplesizes <- c(9999,
                 31740, 37120, 29121, 51396
                 20515, 20623, 17724, 29540,
                 9000, 54986, , 29346,
                 4476, 28130, 16782, 9997)
  


for (i in seq_along(outcomes)) {
  run_mr_pipeline(outcomes[i], ADHD, outcome_sample_size = samplesizes[i])
}




