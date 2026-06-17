# Calculating pairwise variant LD
cd /path/to/your/project
for pheno in $pheno_list; do
    outputpath="./fine_mapping/input_file/ldfile/${pheno}"
    mkdir -p "$outputpath"
    rm $outputpath/*
    while read line; do
        region=`echo $line | awk '{print \$1}'`
        chr=$(echo $region | awk -F':' '{print \$1}')
        snpfile="./fine_mapping/input_file/snplist/${pheno}/${pheno}_${region}.rsid"
        ./software/plink \
            --bfile /path/to/reference/data/Q0_unre_Caucasian_c${chr} \
            --keep ./data/sample_ids.txt \
            --extract ${snpfile} \
            --r inter-chr \
            --ld-window-r2 0 \
            --out ${outputpath}/${pheno}_${region}
    done < ./fine_mapping/pheno_region/${pheno}_merged_regions.txt
    echo "LD file creation Done: ${pheno}"
done

# Running FINEMAP
cd /path/to/your/project/fine_mapping
    /path/to/software/finemap_v1.4.2_x86_64 \
        --sss \
        --log \
        --corr-config 0.95 \
        --n-causal-snps 10 \
        --n-iter 100 \
        --prob-cred-set 0.99 \
        --in-files ./input_file/masterfile/${pheno}_master \
        --n-threads 64
