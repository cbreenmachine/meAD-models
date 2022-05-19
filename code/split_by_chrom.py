import os
import pandas as pd
import argparse

# How to incorporate this script into GNU parallel;
# Need to pass the arguments through a bash script;
# Need to check what's been processed already ala extract script...
# Need to pass input file names (or a directory more likely) but let a bash script handle this

CHROMOSOMES = ['chr' + str(x) for x in range(1,23)]
CHROMOSOMES.append('chrX')
CHROMOSOMES.append('chrY')


def read_and_filter(file_path):
    sample = os.path.basename(file_path).replace('.bed', '')
    df = pd.read_csv(file_path, sep="\t", usecols= ['chr', 'pos', 'methylated', 'unmethylated'], encoding="utf-8")

    df['sample'] = sample
    return df

def split_and_write(df, odir):
    '''Given a dataframe'''
    for c in CHROMOSOMES:
        outfile = os.path.join(odir, c + '.bed')
        
        try:
            if not os.path.isfile(outfile):
                df[df['chr'] == c].to_csv(outfile, sep = '\t', header = 'column_names', index = False)
            else: 
                # else it exists so append without writing the header
                df[df['chr'] == c].to_csv(outfile, sep = '\t', mode = 'a', header = False, index = False)
        except UnicodeDecodeError:
            print(df['sample'][0] + " didn't write")


def main(args):

    if not os.path.exists(odir):
        os.makedirs(odir)

    all_files = []

    for root, subdir, files in os.walk(ROOT):
        if root.endswith("04-extract"):
            for ff in files:
                if ff.endswith(".bed"):
                    all_files.append(os.path.join(root, ff))
    print(all_files)
    print("Working on " + str(len(all_files)) + " files")
    
    counter = 1
    for ff in all_files:
        print("Reading in " + ff)
        df = read_and_filter(ff)
        split_and_write(df)
        print("Finished {0} of {1}".format(counter, len(all_files)))
        counter += 1



if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--idir", default= "../data/04-extract/", help="input file; should be tab separated")
    parser.add_argument("--ofile", default= "../data/05-methylated-coverage-by-chrom/", help = "output file, also tab separated (like BED)")
    args = parser.parse_args()
    main(args)