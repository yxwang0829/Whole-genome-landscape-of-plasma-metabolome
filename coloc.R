library(data.table)
library("coloc")
library(dplyr)

p1 <- 1e-4
p2 <- 1e-4
p12 <- 1e-5

run_coloc_analysis <- function(pheno_snvs_subset_clean, outcome_snvs_subset_clean, r, p1, p2, p12) { 
  coloc_result <- coloc.abf(dataset1=list(
    snp=pheno_snvs_subset_clean$ID,
    beta = pheno_snvs_subset_clean$Meta_Effect,
    varbeta = pheno_snvs_subset_clean$Meta_StdErr^2,
    p=pheno_snvs_subset_clean$Meta_P, 
    MAF=pheno_snvs_subset_clean$Discovery_A1FREQ,
    N= 191681,
    type="quant" ),
    dataset2=list(
      snp=outcome_snvs_subset_clean$id,
      beta=outcome_snvs_subset_clean$beta,
      varbeta=outcome_snvs_subset_clean$se^2,
      p=outcome_snvs_subset_clean$p, 
      MAF=outcome_snvs_subset_clean$maf,
      N= as.numeric(N_finngen),
      type="cc" ),
    p1 = p1, p2 = p2, p12 = p12)
  return(coloc_result)
}

process_coloc_results <- function(coloc_results, metabolite, disease) {
  all_regions_results <- data.frame()
  for (region_name in names(coloc_results)) {
    df_results <- coloc_results[[region_name]]$results
    df_summary <- coloc_results[[region_name]]$summary
    result <- data.frame(nsnps=df_summary['nsnps'],PPH4=df_summary['PP.H4.abf'],
                         coloc_snp = df_results$snp, SNP_PPH4 = df_results$SNP.PP.H4) 
    result <- result %>% slice_max(order_by = SNP_PPH4, n = 1)
    result <- subset(result,result$PPH4>=0.5)
    all_regions_results <- rbind(all_regions_results, result)
  }
  return(all_regions_results)
}

for (i in seq_along(metabolite_list)) {
  m <- metabolite_list[i]
  pheno_gwas <- fread(file_path_m)
  all_metabolite_results <- data.frame()
  for (j in seq_along(disease_list)) {
    d <- disease_list[j]
    outcome_gwas <- fread(file_path_d)
    coloc_results <- list()
    for (r in region_list) {
      region_chromosome <-  d_subset[d_subset$region_identifier==r,]$chr
      region_start <-  d_subset[d_subset$region_identifier==r,]$start_pos
      region_end <-  d_subset[d_subset$region_identifier==r,]$end_pos
      pheno_sub <- pheno_gwas[pheno_gwas$CHR == region_chromosome & 
                                (pheno_gwas$POS >= region_start) & 
                                (pheno_gwas$POS <= region_end)]
      outcome_sub <- outcome_gwas[outcome_gwas$chr == region_chromosome & 
                                    (outcome_gwas$pos >= region_start) & 
                                    (outcome_gwas$pos <= region_end)] 
      common_snps <- intersect(pheno_sub$ID, outcome_sub$id)
      pheno_snvs_subset_clean <- pheno_sub %>%
        distinct(ID, .keep_all = TRUE) %>%
        filter(ID %in% common_snps) %>%
        arrange(ID)
      outcome_snvs_subset_clean <- outcome_sub %>%
        distinct(id, .keep_all = TRUE) %>%
        filter(id %in% common_snps) %>%
        arrange(id)
      coloc_result <- run_coloc_analysis(pheno_snvs_subset_clean, outcome_snvs_subset_clean, r, p1, p2, p12)
      coloc_results[[r]] <- coloc_result 
    }
    disease_results <- process_coloc_results(coloc_results, m, d)
    all_metabolite_results <- rbind(all_metabolite_results, disease_results)
  }
    output_file <- file.path(output_dir, paste0(m, ".csv"))
    write.csv(all_metabolite_results, output_file, row.names = FALSE)
}
