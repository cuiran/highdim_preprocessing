import h5py
import argparse
import pandas as pd
import numpy as np
import pdb

def concat_chr(file_prefix,outfile):
    # concat all chromosomes of ldscores to one large h5 file
    files = [file_prefix+str(i)+'.l2.ldscore.gz' for i in range(1,23)]
    dfs = [pd.read_csv(f,delim_whitespace=True) for f in files]
    df = pd.concat(dfs)
    df.to_csv(outfile+'.csv',sep='\t',index=False)
    with h5py.File(outfile,'w') as f:
        f.create_dataset('dataset',data=df.iloc[:,3:])
    f.close()
    return

def concat_single_chr(ld_list,chrom,outfile):
    lddf = pd.read_csv(ld_list,delim_whitespace=True)
    li = list() # list of dfs
    firstdf = pd.read_csv(lddf.ix[0,'Dir']+chrom+'.l2.ldscore.gz',delim_whitespace=True)
    if lddf.ix[0,'Num_annots']==1:
        firstdf.rename(columns={'L2':lddf.ix[0,'Name']},inplace=True)
    li.append(firstdf)
    mli = list() # list of M files
    mli.append(pd.read_csv(lddf.ix[0,'Dir']+chrom+'.l2.M',delim_whitespace=True,header=None))
    m550li = list() # list of M_5_50 files
    m550li.append(pd.read_csv(lddf.ix[0,'Dir']+chrom+'.l2.M_5_50',delim_whitespace=True,header=None))
    for i in range(1,len(lddf)):
        lddir = lddf.ix[i,'Dir']
        num = lddf.ix[i,'Num_annots']
        name = lddf.ix[i,'Name']
        print('Reading in ldscores for '+name)
        df = pd.read_csv(lddir+chrom+'.l2.ldscore.gz',delim_whitespace=True)
        if num == 1:
            df.rename(columns={'L2':name},inplace=True)
        li.append(df.iloc[:,3:])
        print('Reading in M file for '+name)
        mli.append(pd.read_csv(lddir+chrom+'.l2.M',delim_whitespace=True,header=None))
        print('Reading in M_5_50 file for '+name)
        m550li.append(pd.read_csv(lddir+chrom+'.l2.M_5_50',delim_whitespace=True,header=None))
    print('Concatenating {} LD scores'.format(len(li)))
    allld = pd.concat(li,axis=1)
    allld.to_csv(outfile+'.'+chrom+'.l2.ldscore.gz',sep='\t',index=False,compression='gzip')
    allm = pd.concat(mli,axis=1)
    allm.to_csv(outfile+'.'+chrom+'.l2.M',sep='\t',index=False,header=False)
    allm550 = pd.concat(m550li,axis=1)
    allm550.to_csv(outfile+'.'+chrom+'.l2.M_5_50',sep='\t',index=False,header=False)
    return 

def get_rs(fprefix,fsuffix):
    df = concat_chrs(fprefix,fsuffix)
    return df[['CHR','SNP','BP']]

def concat_to_csv(ld_list,outfile):
    # concatenate all ld scores into one big csv file with SNP names
    lddf = pd.read_csv(ld_list,delim_whitespace=True)
    rsdf = get_rs(lddf.ix[0,'Dir'],'.l2.ldscore.gz')
    li = [rsdf]
    for i in range(len(lddf)):
        lddir = lddf.ix[i,'Dir']
        ldid = lddf.ix[i,'ID']
        print('read in ld score for '+ldid)
        numld = lddf.ix[i,'Num_annots']
        df = concatld_chrs(lddir,ldid,numld)
        li.append(df)
    allld = pd.concat(li,axis=1)
    allld.to_csv(outfile,sep='\t',index=False)
    return

def concatld_chrs(lddir,ldid,numld):
    # given chromosome separated LD files, concatenate them, use ldid as column if numld==1, if numld>1, use the existing column names
    df = concat_chrs(lddir,'.l2.ldscore.gz')
    df.drop(['CHR','SNP','BP'],axis=1,inplace=True)
    if numld==1:
        df.columns=[ldid]
    return df

def concat_chrs(fprefix,fsuffix):
    files = [fprefix + str(i) + fsuffix for i in range(1,23)]
    dfs = [pd.read_csv(f,delim_whitespace=True) for f in files]
    df = pd.concat(dfs)
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--concat_chr',action='store_true')
    parser.add_argument('--ldfiles')
    parser.add_argument('--outfile',help="output prefix not including dot")
    parser.add_argument('--chrom')
    parser.add_argument('--ld_list')
    parser.add_argument('--combine_onechr',action='store_true')
    parser.add_argument('--to_csv',action='store_true')
    parser.add_argument('--rs_ref')
    args = parser.parse_args()

    if args.concat_chr:
        concat_chr(args.ldfiles,args.outfile)
    elif args.combine_onechr:
        concat_single_chr(args.ld_list,args.chrom,args.outfile)
    elif args.to_csv:
        concat_to_csv(args.ld_list,args.outfile)
