#!/bin/bash

input_file="/path/to/input/file/total_validated_meta_sig_metwgs_single_variant.tsv"
output_dir="/path/to/output/directory/clump"

mkdir -p "${output_dir}"
metabolites=$(tail -n +2 "${input_file}" | cut -f1 | sort | uniq)

for metabolite in ${metabolites}; do
    echo "Processing metabolite: ${metabolite}"

    temp_input="${output_dir}/${metabolite}_temp.tsv"
    head -n 1 "${input_file}" > "${temp_input}"
    grep "^${metabolite}[[:space:]]" "${input_file}" >> "${temp_input}"

    metabolite_dir="${output_dir}/${metabolite}"
    mkdir -p "${metabolite_dir}"

    for chr in {1..22}; do
        echo "Processing chromosome ${chr} for ${metabolite}"
        /path/to/software/plink2 \
        --bfile /path/to/genetic/data/c${chr} \
        --rm-dup force-first \
        --keep /path/to/keep_file.txt \
        --clump "${temp_input}" \
        --clump-p1 1.6e-10 \
        --clump-p2 1.6e-10 \
        --clump-r2 0.1 \
        --clump-kb 1000 \
        --clump-field Meta_P-value \
        --clump-snp-field ID \
        --chr ${chr} \
        --out "${metabolite_dir}/clumping_results_chr${chr}"
    done
    rm "${temp_input}"
done

echo "Clumping process completed for all metabolites and chromosomes."
