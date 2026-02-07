library(dplyr)
library(data.table)
library("TwoSampleMR")
library("readr")

for (j in seq_along(disease_list)) {
  d <- disease_list[j]
  outcome_data<-read_outcome_data(
    filename = file_path_d,    
    sep = "\t",
    snp_col = "id",
    effect_allele_col = "EA", 
    other_allele_col = "NEA",
    pval_col = "p",
    chr_col = "chr",
    beta_col = "beta", 
    se_col = "se",
    eaf = "maf")
  
  for (i in seq_along(metabolite_list)) {
    m <- metabolite_list[i]
    exposure_data<-read_exposure_data(
      filename = leadsnp_file,
      sep = "\t",
      snp_col = "ID",        
      beta_col = "beta", 
      se_col = "se", 
      effect_allele_col = "Allele1", 
      other_allele_col = "Allele2", 
      pval_col = "p",
      eaf_col = "A1_FREQ"
    ) 

    dat <- harmonise_data(exposure_dat=exposure_data, outcome_dat=outcome_data) 
    if (nrow(dat) == 0) {
      next
    } else if (nrow(dat) == 1) {
      results <- mr(dat, method_list=c("mr_wald_ratio"))
      results.withOR <- generate_odds_ratios(results)
    } else if (nrow(dat) > 1) {
      results <- mr( dat)
      results.withOR <- generate_odds_ratios(results)
      dat$samplesize.outcome<-No  
      dat$samplesize.exposure<-Ne
      dat.mr_heterogeneity<-mr_heterogeneity( dat)   
      dat.mr_pleiotropy_test<-mr_pleiotropy_test( dat) 
      leaveoneout <- mr_leaveoneout( dat)
      leaveoneout.withOR <- generate_odds_ratios( leaveoneout)
      single_snp_analysis <- mr_singlesnp(dat) 
      results.single_snp_analysis.withOR <- generate_odds_ratios(single_snp_analysis)
    }
  }
}