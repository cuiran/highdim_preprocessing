import pandas as pd
import argparse
from os import path
import pdb

def filter_annots(annot_list,chrom,annot_bim,target_bim,output_prefix):
    '''
    Filter annotations for more SNPs to less SNPs.

    annot_list is a file with extensino .annotcts, with columns 
    'ID','Name','Dir','Num_annots','Thin'

    annot_bim is the bim file that corresponding to the SNPs in the annotations

    target_bim is the bim file that you want to filter the SNPs down to.

    output_prefix specifies the folder you want to store the filtered annotations in
    '''
    # read in annotation bim file if exist
    if annot_bim is not None:
        bimdf = pd.read_csv(annot_bim+chrom+'.bim',delim_whitespace=True,header=None,usecols=[1])
        bimdf.columns = ['SNP']
    # read in target bim file
    bimdf_filt = pd.read_csv(target_bim+chrom+'.bim',usecols=[1],delim_whitespace=True,header=None)
    bimdf_filt.columns=['SNP']
    snps = bimdf_filt['SNP'].tolist()
    print('There are {} SNPs in the filtered annotation'.format(len(snps)))
    # process every annotation in the list
    list_df = pd.read_csv(annot_list,delim_whitespace=True)
    for i in range(len(list_df)):
        name = list_df.loc[i,'Name']
        prefix = list_df.loc[i,'Dir']
        thin = list_df.loc[i,'Thin']
        print('Filtering annotation for cell type {}'.format(name))
        if path.exists(prefix+chrom+'.annot.gz'):
            annot_df = pd.read_csv(prefix+chrom+'.annot.gz',delim_whitespace=True)
            print('There are {} SNPs in the original annotation'.format(annot_df.shape[0]))
            if thin:
                annot_df = pd.concat([bimdf,annot_df],axis=1)
            annot_df.set_index('SNP',inplace=True)
            annot_filt = annot_df.loc[snps,:]
            if not thin:
                annot_filt.reset_index(inplace=True)     
            annot_filt.to_csv(output_prefix+name+'.'+chrom+'.annot.gz',sep='\t',index=False,compression='gzip') 
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--filter-annots',action='store_true')
    parser.add_argument('--chrom',type=str)
    parser.add_argument('--annot-list')
    parser.add_argument('--annot-bim')
    parser.add_argument('--target-bim')
    parser.add_argument('--output-prefix')
    args = parser.parse_args()

    if args.filter_annots:
        filter_annots(args.annot_list,args.chrom,args.annot_bim,args.target_bim,args.output_prefix)
