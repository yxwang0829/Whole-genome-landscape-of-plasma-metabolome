library(data.table)
library(dplyr)
library(combinat)
library(hash)
library(car)
library(ggplot2)
library(ggrepel)

   
for(d in start_index:end_index) {
  outcome_name=unique(index$disease)[d]   
  matrix_beta_filtered <- matrix_beta[matrix_beta$ID %in% unique(lead_subset$ID), c("ID", ..met_list)]
  exposure_matrix <- as.matrix(matrix_beta_filtered[, -1])  
  
  cor_matrix <- cor(exposure_matrix)
  diag(cor_matrix) <- NA
  high_cor_pairs <- which(abs(cor_matrix) > 0.99, arr.ind = TRUE)
  removed_metabolites <- logical(ncol(exposure_matrix))    
  if(nrow(high_cor_pairs) > 0) {
    for (pair in 1:nrow(high_cor_pairs)) {
      i <- high_cor_pairs[pair, 1]
      j <- high_cor_pairs[pair, 2]
      if (removed_metabolites[i] | removed_metabolites[j]) next
      if (!(colnames(exposure_matrix)[i] %in% non_derived_list) && 
          !(colnames(exposure_matrix)[j] %in% non_derived_list)) {
        removed_metabolites[i] <- TRUE
      } else {
        if (!(colnames(exposure_matrix)[i] %in% non_derived_list)) {
          removed_metabolites[i] <- TRUE
        } else {
          removed_metabolites[j] <- TRUE
        }
      }
    }

    if(any(removed_metabolites)) {
      exposure_matrix <- exposure_matrix[, !removed_metabolites]
      met_list <- colnames(exposure_matrix)
      lead_subset <- subset(lead, lead$Metabolite %in% met_list)
      matrix_beta_filtered <- matrix_beta[matrix_beta$ID %in% unique(lead_subset$ID), 
                                          c("ID", ..met_list)]
      exposure_matrix <- as.matrix(matrix_beta_filtered[, -1])
    }
  }
  snps_list <- matrix_beta_filtered$ID
  
  matrix_se_filtered <- matrix_se[matrix_se$ID %in% unique(lead_subset$ID), c("ID", ..met_list)]
  exposure_se_matrix <- as.matrix(matrix_se_filtered[, -1]) 
  
  matrix_metabolite_filtered = matrix_metabolite[, c("rn", ..met_list)] 
  matrix_metabolite_filtered = matrix_metabolite_filtered[match(met_list, matrix_metabolite_filtered$rn),]
  Pcov=as.matrix(matrix_metabolite_filtered[,-1])
  
  outcome_data=left_join(matrix_beta_filtered[,1],outcome)[,c('beta','se','EA','NEA','maf')]
  complete_cases <- complete.cases(outcome_data)
  exposure_matrix_complete <- exposure_matrix[complete_cases,]
  exposure_se_matrix_complete <- exposure_se_matrix[complete_cases,]
  outcome_data <- outcome_data[complete_cases,]
  snps_list <- snps_list[complete_cases]
  
  mr_input <- new("mvMRInput",
                  betaX = exposure_matrix_complete,
                  betaY = as.matrix(outcome_data$beta),
                  betaXse = exposure_se_matrix_complete,     
                  betaYse = as.matrix(outcome_data$se), 
                  exposure = met_list,
                  outcome = outcome_name,
                  snps = snps_list,
                  effect_allele = outcome_data$EA, 
                  other_allele = outcome_data$NEA,  
                  eaf = outcome_data$maf,        
                  correlation = Pcov)
  results <- mr_bma(
    mvmr_input = mr_input,
    prior_prob = 0.1,
    prior_sigma = 0.25,
    kmin=1,
    kmax=min(5, length(met_list)),
    remove_outliers = FALSE,
    remove_influential = FALSE,
    calculate_p = FALSE,
    nrepeat = 10000
  )
  model_best <- results$model_best
  mip_table <- results$mip_table
  mrbma_output <- results$mrbma_output
  outlier_res <- results$outlier_res
  influential_res <- results$influential_res
  bma_results <- sss.report.mr.bma(mrbma_output, 
                                   top = length(met_list), 
                                   write.out = TRUE,
                                   csv.file.name = file.path(outputpath, "bma_results.csv"))
  best_model_results <- sss.report.best.model(mrbma_output, 
                                              top = length(met_list))
  gc()
}

