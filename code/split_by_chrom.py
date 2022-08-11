import os
import pandas as pd
import argparse
import logging

# How to incorporate this script into GNU parallel;
# Need to pass the arguments through a bash script;
# Need to check what's been processed already ala extract script...
# Need to pass input file names (or a directory more likely) but let a bash script handle this

# Constants
CHROMOSOMES = ['chr' + str(x) for x in range(1,23)]
CHROMOSOMES.append('chrX')
CHROMOSOMES.append('chrY')

def make_dir(new_dir):
    '''create the directory if it does not exist'''
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)


def read_and_filter(file_path):
    '''Given a path to bed file (e.g. ../data/04-extract/100.bed) reads in 
    as dataframe and adds column with sample id'''
    # Sample is file name without relative dir and extension
    sample = os.path.basename(file_path).replace('.bed', '')

    # Read in
    df = pd.read_csv(file_path, sep="\t", encoding="utf-8")

    # Add column
    df['sample'] = sample
    return df

def split_and_write(df, odir, counter):
    '''Given a dataframe, writes out by chromosome'''
    # Use counter in case we want to overwrite old files
    for c in CHROMOSOMES:
        outfile = os.path.join(odir, c + '.bed')
        
        if counter == 1:
            df[df['chrom'] == c].to_csv(outfile, sep = '\t', header = 'column_names', index = False)
        else: 
            # else it exists so append without writing the header
            df[df['chrom'] == c].to_csv(outfile, sep = '\t', mode = 'a', header = False, index = False)



def main(args):

    idir, odir = args.idir, args.odir
    make_dir(odir)    

    all_files = []

    # Stupid python directory crap
    for (dirpath, dirnames, filenames) in os.walk(idir):
        for ff in filenames:
            if ff.endswith(".bed"):
                all_files.append(os.path.join(idir, ff))
    print("Working on " + str(len(all_files)) + " files")
    
    counter = 1
    for ff in all_files:
        print("Reading in " + ff)

        df = read_and_filter(ff)
        split_and_write(df, odir, counter)
        
        # Console output
        print("Finished {0} of {1}".format(counter, len(all_files)))
        counter += 1



if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--idir", default= "../data/04-extract/", help="input file; should be tab separated")
    parser.add_argument("--odir", default= "../data/06-counts-by-chrom/", help = "output file, also tab separated (like BED)")
    args = parser.parse_args()
    main(args)