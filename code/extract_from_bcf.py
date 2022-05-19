# extract_from_bcf.py
# Broadly speaking, takes 

# Piped from bcftools query:
# %CHROM\t%POS\t[%CS]\t%REF\t[%MC8{4}]\t[%MC8{5}]\t[%MC8{6}]\t[%MC8{7}]\n
# In other words:
# chrom,pos,strand,reference,num_A,num_C,num_G,num_T

# Takes input with a line that looks like
# chr1    15941   16341   +       7       7      
# chr1    15942   16342   -       8       9      
# And merges across strand. In other words, produces
# chr1    15941   16341   +/-     15      16    

from sys import stdin, stdout
import pandas as pd
import numpy as np
import argparse

# See the bcftools query command in ./process.sh extract
COLNAMES=["chrom","chromStart","strand","reference","num_A","num_C","num_G","num_T"] 



def read_data(path, col_names=COLNAMES):
    '''read as pandas data frame and assign column names
    
    '''
    
    df = pd.read_csv(path, sep='\t', names=col_names, header=None)
    df['chromStart'] -= 1 # BCF to BED coordinate system

    # Cleaning
    df.dropna(inplace = True)
    df.drop_duplicates(inplace = True)
    df['chromEnd'] = df['chromStart'] + 2
    return(df)



def merge_strands(data):
    '''Merges (un)methylated reads across strands
    '''
    df = data.copy() # stupid python nonsense

    # Move the index of negative stranded CpGs back one
    df.loc[df['strand'] == '-', ['chromStart', 'chromEnd']] -= 1
    # Then group by genomic coordinates, and merge rows 
    # methylated and coverage get summed
    # strand is turned into +/-
    # And sequence comes from first row in grouping (top strand CpG)
    df = (df.groupby(by = ['chrom', 'chromStart', 'chromEnd']) \
    .agg({'strand': '/'.join, 'methylated': 'sum', 'unmethylated': 'sum', 'coverage': 'sum'})) #, 'sequence': lambda x: x.iloc[0]}))

    return(df)


def main(args):
    df = read_data(args.ifile)

    df['methylated'] = np.where(df['strand'] == "+", df['num_C'], df['num_G'])
    df['unmethylated'] = np.where(df['strand'] == "+", df['num_T'], df['num_A'])
    df['coverage'] = df['methylated'] + df['unmethylated']
    

    if args.merge:
        df = merge_strands(df)

    df.drop(["unmethylated"], axis=1, inplace=True)
    df.to_csv(args.ofile, sep='\t')



if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ifile", default= stdin, help="input file; should be tab separated")
    parser.add_argument("--ofile", default= stdout, help = "output file, also tab separated (like BED)")
    parser.add_argument("--merge", action="store_true", help = "")
    args = parser.parse_args()
    main(args)
    