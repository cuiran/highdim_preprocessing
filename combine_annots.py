import h5py
import argparse
import pandas as pd
import numpy as np
import pdb

def concat_chr(file_prefix,outfile):
    # concat all chromosomes of ldscores to one large h5 file
    files = [file_prefix+str(i)+'.annot.gz' for i in range(1,23)]
    dfs = [pd.read_csv(f,delim_whitespace=True) for f in files]
    df = pd.concat(dfs)
    df.to_csv(outfile+'.csv',sep='\t',index=False)
    with h5py.File(outfile,'w') as f:
        f.create_dataset('dataset',data=df.iloc[:,4:])
    f.close()
    return

def concat_single_chr(ld_list,chrom,outfile):
    lddf = pd.read_csv(ld_list,delim_whitespace=True)
    li = list() # list of dfs
    firstdf = pd.read_csv(lddf.ix[0,'Dir']+chrom+'.annot',delim_whitespace=True)
    if lddf.ix[0,'Num_annots']==1: 
        firstdf.rename(columns={'ANNOT':lddf.ix[0,'Name']},inplace=True)
    li.append(firstdf)
    for i in range(1,len(lddf)):
        lddir = lddf.ix[i,'Dir']
        num = lddf.ix[i,'Num_annots']
        name = lddf.ix[i,'Name']
        print('Reading in annotations for '+name)
        df = pd.read_csv(lddir+chrom+'.annot',delim_whitespace=True)
        if num==1:
            df.rename(columns={'ANNOT':name},inplace=True)
        if lddf.ix[i,'Thin'] == True:
            li.append(df.iloc[:,:])
        elif lddf.ix[i,'Thin'] == False:
            li.append(df.iloc[:,4:])
    allann = pd.concat(li,axis=1)
    print(allann.shape)
    allann.to_csv(outfile+'.'+chrom+'.annot.gz',index=False,sep='\t',compression='gzip')
    return 

def corr(args):
    annot_prefix = args.annot_prefix
    outfile = args.outfile

    anns = list()
    for i in range(1,23):
        chrom = str(i)
        print('chromosome '+chrom)
        annfname = annot_prefix+chrom+'.annot.gz'
        df = pd.read_csv(annfname,delim_whitespace=True,header=None)
        if all(x in df.columns for x in ['CHR','BP','SNP','CM']):
            df.drop(['CHR','BP','SNP','CM'],axis=1,inplace=True)
        d = np.array(df.iloc[:,:])
        anns.append(d)
    ann_matrix = np.concatenate(anns,axis=0)
    corr_matrix = np.corrcoef(ann_matrix,rowvar=False)
    corrdf = pd.DataFrame(corr_matrix,columns=range(corr_matrix.shape[1]))
    corrdf.to_csv(outfile,sep='\t',index=False,header=False)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--concat_chr',action='store_true')
    parser.add_argument('--annotfiles')
    parser.add_argument('--outfile')
    parser.add_argument('--chrom')
    parser.add_argument('--annot_list')
    parser.add_argument('--combine_onechr',action='store_true')
    parser.add_argument('--corr',action='store_true')
    parser.add_argument('--annot_prefix')
    args = parser.parse_args()

    if args.concat_chr:
        concat_chr(args.annotfiles,args.outfile)
    elif args.combine_onechr:
        concat_single_chr(args.annot_list,args.chrom,args.outfile)
    elif args.corr:
        corr(args)
