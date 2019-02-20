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
    for i in range(len(lddf)):
        lddir = lddf.ix[i,'Dir']
        df = pd.read_csv(lddir+chrom+'.annot.gz',delim_whitespace=True)
        if lddf.ix[i,'Thin'] == 'T':
            li.append(df.iloc[:,:])
        elif lddf.ix[i,'Thin'] == 'F':
            li.append(df.iloc[:,4:])
        print(lddir)
        print(df.iloc[:,4:].shape)
    allann = np.concatenate(li,axis=1)
    print(allann.shape)
    #f = h5py.File(outfile,'w')
    #f.create_dataset('dataset',data=allld)
    #f.close()
    df = pd.DataFrame(data=allann,columns=range(allann.shape[1]))
    df.to_csv(outfile,header=False,index=False,sep='\t',compression='gzip')
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