##REGENIE
#https://rgcgithub.github.io/regenie/options
regenie \
  --step 1 \
  --bed /path/to/your/data/bed/file \
  --extract /path/to/your/data/snplist \
  --keep /path/to/your/data/${set}_sample_ids.txt \
  --phenoFile /path/to/your/data/phenofile.txt \
  --covarFile /path/to/your/data/covariatefile.txt \
  --catCovarList Sex,Ethnicity,Month,Center,Provider,Batch \
  --maxCatLevels 26 \
  --apply-rint \
  --force-qt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix /path/to/your/temp/preds \
  --out /path/to/your/output/${set}_step1_output
