#!/bin/bash

python /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/pyscripts/combine_ld.py \
    --to_csv \
    --ld_list $1 \
    --outfile /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annotations/EAS/${2}/${3}_combined_ld.csv
