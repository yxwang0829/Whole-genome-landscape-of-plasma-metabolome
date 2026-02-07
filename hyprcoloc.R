library(data.table)
library(hyprcoloc)
library(dplyr)
library(stringr)

all_results <- data.frame()
for(chr in names(chr_groups)) {  
  current_index <- chr_groups[[chr]]
  for (i in seq_len(nrow(current_index))) { 
    current_row <- current_index[i,]  
      for(t in tissues){ 
        eqtl_data <- eqtl_data_list[[t]] %>% filter(Probe == current_row$probe)
        combined_data <- metabolite_data %>%
          inner_join(disease_data, by = "SNP") %>%
          inner_join(eqtl_data[, c('SNP', 'beta_e', 'se_e')], by = "SNP") %>%
          na.omit()
        
        betas <- as.matrix(combined_data[, c('beta_m', 'beta_d', 'beta_e')])
        ses <- as.matrix(combined_data[, c('se_m', 'se_d', 'se_e')])
        traits <- c(current_row$metabolite, current_row$disease, current_row$probe)
        rsid <- combined_data$SNP
        
        result <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid, binary.outcomes = c(0,1,0))
        all_results <- rbind(all_results, result_df)
        }
      } 
  } 


