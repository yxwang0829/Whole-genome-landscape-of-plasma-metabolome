## Calculate inflation factor lambda
rm(list=ls())
suppressMessages(library(data.table))
library(qqman)
setwd("/path/to/your/project")  # Set your project directory
phenolist <- names(as.data.frame(fread("./data/phenofile.txt", header = T)))[-c(1, 2)]
result <- data.table(Metabolite = phenolist,
                     Lambda = 0)
for (i in 1:length(phenolist)) {
  pheno <- phenolist[i]
  path <- paste0("./meta/single_summary/sum_", pheno, "_meta_single_variant_wgs_summary_statistics.tsv")
  p <- as.numeric(unlist(fread(path, select = "Meta_P-value", header = T)))
  # Inflation factor
  z <- qnorm(p/2)
  result[i, "Lambda"] <- round(median(z^2, na.rm = TRUE) / 0.454, 3)
  print(i)
}
writepath <- "./meta/single_lambda.txt"
fwrite(result, writepath, sep = "\t", quote = F, row.names=F, col.names = T)
