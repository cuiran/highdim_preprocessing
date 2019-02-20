#!/bin/bash

# combine all annots into chromosome separated annot.gz files
for i in {1..22}
do sbatch -p short -t 0-12:00 --mem=8000 -o ../output/combine_annot_chr${i}.out -e ../output/combine_annot_chr${i}.err \
combine_annots.sh $i /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/ref_files/allannot_dirlist_eas.csv all715 all
done

# combine all ldscores into one csv file
sbatch -p short -t 0-12:00 --mem=32000 -o ../output/combine_ldscore.out -e ../output/combine_ldscore.err \
combine_ldscore.sh /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/ref_files/allannot_dirlist_eas.csv all715 all

# combine all baselineLD annots
for i in {1..22}
do sbatch -p short -t 0-12:00 --mem=8000 -o ../output/combine_blannot_chr${i}.out -e ../output/combine_blannot_chr${i}.err \
combine_annots.sh $i /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/ref_files/baselineLD_dirlist_eas.csv baselineLD_v2.1 BL
done 

# combine all baselineLD ldscores into one csv file
sbatch -p short -t 0-12:00 --mem=32000 -o ../output/combine_ldscore.out -e ../output/combine_blldscore.err \
combine_ldscore.sh /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/ref_files/baselineLD_dirlist_eas.csv baselineLD_v2.1 BL
