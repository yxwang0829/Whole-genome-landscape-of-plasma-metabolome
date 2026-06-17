rm(list = ls())
library(data.table)
library(dplyr)
library(mediation)
setwd('/path/to/project/directory')
sig_phewas_asso <- fread("./path/to/sig_ketone_leadsnp_phewas_results.txt", header = T)
lead_asso <- fread("./path/to/ketone_bodies_single_anno_pool_together_adj_glucose.txt", select = c('ID', 'Metabolite'), header = T)
sig_phewas_asso <- sig_phewas_asso %>%
  as.data.frame() %>%
  distinct(ID, endpoint) %>%
  left_join(lead_asso, by = 'ID', relationship = "many-to-many")
metabolome_data <- fread('./data/Caucasian_regenie_phenofile.txt', header = T) %>% as.data.frame() %>% dplyr::select(-FID)
covar_file <- fread('./data/Caucasian_regenie_corvarfile.txt', header = T) %>% dplyr::select(-FID)
covariates <- names(covar_file)[-c(1, 2)]
info <- c("ACME", "ACME_P", "ADE", "ADE_P", "Total_Effect", "Total_Effect_P", 
          "N", "Prop_mediated", "Prop_mediated_P", "ACME_LCI", "ACME_UCI", 
          "ADE_LCI", "ADE_UCI", "Total_Effect_LCI", "Total_Effect_UCI", 
          "Prop_mediated_LCI", "Prop_mediated_UCI")
for (i in 1:nrow(sig_phewas_asso)) {
  variant <- sig_phewas_asso[i, 'ID']
  chr <- as.numeric(sub(":.*", "", variant))
  metabolite <- sig_phewas_asso[i, 'Metabolite']
  endpoint <- sig_phewas_asso[i, 'endpoint']
  snp_raw_file <- as.data.frame(fread(paste0('./path/to/snp_exp/chr', chr, '_single_snp_mutation_info.raw'), header = T))
  colnames(snp_raw_file) <- lapply(colnames(snp_raw_file), function(x) gsub("_.*$", "", x))
  snp_raw_file <- snp_raw_file %>% dplyr::select(IID, exposure = all_of(variant)) %>% na.omit()
  endpoint_file <- fread('./path/to/endpoint_file.FullDisease.txt', select = c('IID', endpoint), header = T) %>%
    filter(IID %in% metabolome_data$IID) %>%
    na.omit()
  path <- paste0('./path/to/endpoint/, 'step2_Caucasian_chr', chr, '_', endpoint, '.regenie.ids')
  endpoint_id <- fread(path, header = T)$IID
  data <- snp_raw_file %>%
    inner_join(metabolome_data %>% dplyr::select(IID, all_of(metabolite)) %>% na.omit(), by = 'IID') %>%
    inner_join(covar_file, by = 'IID') %>%
    inner_join(endpoint_file, by = 'IID') %>%
    filter(IID %in% endpoint_id)
  print(paste('i: ', i, ';  exposure: ', variant, ', mediator: ', metabolite, ', endpoint: ', endpoint))
  tryCatch(
    {
      formula_mediator <- paste0(metabolite, " ~ exposure + ", paste(covariates, collapse = " + ")) %>% as.formula()
      formula_outcome <- paste0(endpoint, " ~ exposure + ", metabolite, " + ", paste(covariates, collapse = " + ")) %>% as.formula()
      
      model_mediator <- lm(formula_mediator, data = data)
      model_outcome <- glm(formula_outcome, data = data,
                           family = binomial(link="probit"))
      set.seed(123)
      mediation_result <- mediate(
        model_mediator, model_outcome,
        treat = 'exposure', mediator = metabolite, outcome = endpoint,
        boot = TRUE, sims = 1000
      )
      mediation_summary <- mediation_result %>% summary()
      
      sig_phewas_asso[i, info] <- c(
        mediation_summary$d.avg, mediation_summary$d.avg.p, mediation_summary$z.avg, mediation_summary$z.avg.p,
        mediation_summary$tau.coef, mediation_summary$tau.p, mediation_summary$nobs, mediation_summary$n.avg, mediation_summary$n.avg.p,
        mediation_summary$d.avg.ci[1], mediation_summary$d.avg.ci[2], mediation_summary$z.avg.ci[1], mediation_summary$z.avg.ci[2],
        mediation_summary$tau.ci[1], mediation_summary$tau.ci[2],
        mediation_summary$n.avg.ci[1], mediation_summary$n.avg.ci[2]
      )
    },
    error = function(e) {
      print(paste("Error in iteration ", i, ": ", e$message))
    }
  )
}
writepath <- "./ketone_bodies/sig_ketone_phewas_asso_mediation_result.txt"
fwrite(sig_phewas_asso, writepath, sep = "\t", row.names = F, col.names = T, quote = F)
