#!/bin/bash

python /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/pyscripts/combine_annots.py \
    --combine_onechr --annot_list $2 --chrom $1 --outfile ../annotations/EAS/${3}/${4}.${1}.annot.gz
