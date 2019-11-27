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

    ALLELE0 is the reference allele, while ALLELE1 is the effect allele.
    '''
    # use munge_sumstats to format sumstats
    subprocess.call(['python',path_to_ldsc+'munge_sumstats.py',
    '--sumstats',ss_fname,
    '--N',N,
    '--out',output_prefix,
    '--p','P_LINREG',
    '--a1','ALLELE1',
    '--a2','ALLELE0'])
    return

def formatSAIGEss(ss_fname,N,output_prefix,path_to_ldsc):
    ''' 
    Format SAIGE summary statistics files.

    Such file should have columns (in order): 'v', 'CHR', 'POS', 'SNPID', 'Allele1', 'Allele2', 'AC_Allele2',
    'AF_Allele2', 'N', 'BETA', 'SE', 'Tstat', 'p.value', 'p.value.NA',
    'Is.SPA.converge', 'varT', 'varTstar', 'N.Cases', 'N.Controls',
    'AF.Cases', 'AF.Controls'

    Use p.value as the --p flag input

    Allele1 is the reference allele, while Allele2 is the effect allele
    '''
    subprocess.call(['python',path_to_ldsc+'munge_sumstats.py',
    '--sumstats',ss_fname,
    '--N',N,
    '--out',output_prefix,
    '--p','p.value',
    '--a1','Allele2',
    '--a2','Allele1',
    '--snp','SNPID'])
    return

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ss-fname')
    parser.add_argument('--N')
    parser.add_argument('--output-prefix')
    parser.add_argument('--path-to-ldsc')
    parser.add_argument('--bolt',action='store_true')
    parser.add_argument('--saige',action='store_true')
    args = parser.parse_args()

    if args.bolt:
        formatBOLTss(args.ss_fname,args.N,args.output_prefix,args.path_to_ldsc)
    elif args.saige:
        formatSAIGEss(args.ss_fname,args.N,args.output_prefix,args.path_to_ldsc) 
