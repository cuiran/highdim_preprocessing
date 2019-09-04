import pandas as pd
import argparse
import pdb
import subprocess

def formatBOLTss(ss_fname,N,output_prefix,path_to_ldsc):
    '''
    Format BOLT summary statistics files.

    Such file should have columns (in order): 'v', 'SNP', 'CHR', 'BP', 'GENPOS', 'ALLELE1', 'ALLELE0', 'A1FREQ',
    'INFO', 'CHISQ_LINREG', 'P_LINREG', 'BETA', 'SE', 'CHISQ_BOLT_LMM_INF',
    'P_BOLT_LMM_INF'

    Use 'CHISQ_LINREG' as the chisq statistics

    Resulting formatted summary statistics should have ['SNP','A1','A2','Z','N']
    '''
    # remove SNP column
    temp_fname = ss_fname+'.temp'
    f = open(temp_fname,"w")
    subprocess.call(['cut','-f1-1,3-',ss_fname],stdout=f)
    f.close()
    # use munge_sumstats to format sumstats
    subprocess.call(['python',path_to_ldsc+'munge_sumstats.py',
    '--sumstats',temp_fname,
    '--N',N,
    '--out',output_prefix,
    '--p','P_LINREG',
    '--a1','ALLELE1',
    '--a2','ALLELE0',
    '--snp','v'])
    return

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ss-fname')
    parser.add_argument('--N')
    parser.add_argument('--output-prefix')
    parser.add_argument('--path-to-ldsc')
    parser.add_argument('--bolt',action='store_true')
    args = parser.parse_args()

    if args.bolt:
        formatBOLTss(args.ss_fname,args.N,args.output_prefix,args.path_to_ldsc)

    
