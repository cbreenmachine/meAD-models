import pandas as pd
import argparse
import os
import numpy as np
from sklearn.decomposition import IncrementalPCA

def read_and_join_all(all_files):
    df_raw = pd.concat((pd.read_csv(f, sep="\t") for f in all_files), ignore_index=True)
    return df_raw

def impute_col_mean(df, columns):
    '''Handles
    
    Inputs:
    df : pandas dataframe
        A data frame (usually coverage or methylated reads)
    
    No smoothing because that's handled in DSS later on.
    '''
    # Adapted from https://stackoverflow.com/questions/33058590/pandas-dataframe-replacing-nan-with-row-average
    m = df[columns].mean(axis=1, skipna=True).round()
    for c in columns:
        df[c].fillna(m, inplace=True)


def run_pca(X, n_components = 100):
    '''computes principal components incrementally
    '''
    # Makes the computation take a reasonable amount of time
    ipca = IncrementalPCA(n_components = n_components, batch_size=10000)
    X_ipca = ipca.fit_transform(X)
    
    pcs = pd.DataFrame(X_ipca)
    pcs.columns = ['PC' + str(x+1) for x in range(n_components)]

    data = {'var_explained': ipca.explained_variance_ratio_, 
            'PCs' : [x+1 for x in range(n_components)]}
    var = pd.DataFrame(data)
    return pcs, var


def make_name(idir, chrom, sample_num):
    '''inDirectory/chr18.100.bed'''
    return os.path.join(idir, chrom + "." + str(sample_num) + ".bed")


def main(args):

    # Two seperate folders
    if args.drop_MCI:
        odir = os.path.join(args.odir, "withoutMCI")
    else:
        odir = os.path.join(args.odir, "withMCI")

    # Handle directory structure
    if not os.path.exists(odir):
        os.makedirs(odir)

    idir, samples_file, chrom = args.idir, args.samples_file, args.chrom
    samples_df = pd.read_csv(samples_file)

    valid_files = []
    for x in samples_df['sample_id']:
        ff = make_name(idir, chrom, x)
        if os.path.exists(ff):
            valid_files.append(ff)
    
    # Simplifies the loading/PC computation to speed up
    if args.dry_run:
        all_files = valid_files[0:10]
    else:
        all_files = valid_files

    # This takes awhile
    df_raw = read_and_join_all(all_files).drop_duplicates(['sample', 'chromStart'])
    df = df_raw.pivot(index='chromStart', columns='sample', values=['methylated', 'coverage'])
    
    # Split into coverage and methylated dataframes 
    # (Probably a way to do this with multi-indexing, but this is easier to understand
    # and more R-like)
    # reset_index gives chromStart as a column again, so hold off until imputation is done
    Cov = df['coverage'].copy(deep=True).rename_axis(None, axis=1)
    M = df['methylated'].copy(deep=True).rename_axis(None, axis=1)
    
    Cov.columns = Cov.columns.astype(str)
    M.columns = M.columns.astype(str)

    loaded_samples = set(M.columns.to_list())
    
    # DO NOT CHANGE; the underlying optimization problem 
    # breaks down when you approach 100% variance explained 
    # I.E. number of components == number of samples breaks
    n_components = int(len(loaded_samples) / 4)

    print("Imputing missing values...")
    for group in list(set(samples_df['diagnostic_group'])):
        all_group_samples = set([str(x) for x in samples_df[samples_df.diagnostic_group == group]['sample_id']])
        cols = list(all_group_samples.intersection(loaded_samples))
        
        if args.drop_MCI and group == "MCI":
            continue
        else:
            print("Group {0} has {1} samples".format(group, len(cols)))
            impute_col_mean(M, cols)
            impute_col_mean(Cov, cols)
      
    print("Determining drop positions")
    ix_1 = np.where(Cov.isna().sum(axis=1) > 0)[0]
    ix_2 = np.where(M.isna().sum(axis=1) > 0)[0]
    ix_3 = np.where((Cov == 0).sum(axis=1) > 0)[0]

    print("Resetting indices")
    Cov.reset_index(inplace=True)
    M.reset_index(inplace=True)  

    ix_to_drop = list(set(np.concatenate([ix_1, ix_2, ix_3])))
    print(str(len(ix_to_drop)) + " of " + str(len(Cov)) + " positions will be dropped. ")

    print("Dropping")
    Cov.drop(ix_to_drop, axis=0, inplace=True)
    M.drop(ix_to_drop, axis=0, inplace=True)

    ofile = os.path.join(odir, chrom + ".M.bed")
    M.to_csv(ofile, sep='\t', index=False)

    ofile = os.path.join(odir, chrom + ".Cov.bed")
    Cov.to_csv(ofile, sep='\t', index=False)

    # Finally, run PCA
    X = M.drop('chromStart', axis=1).div(Cov.drop('chromStart', axis=1)).transpose()
    X_centered = X - X.mean(axis=0)
    pcs, var = run_pca(X_centered, n_components)

    pcs['sample_id'] = Cov.columns.to_list()[1:]

    ofile = os.path.join(odir, chrom + ".PCs.csv")
    pcs.to_csv(ofile, index=False)

    ofile = os.path.join(odir, chrom + ".PCs.varexp.csv")
    var.to_csv(ofile, index=False)

  
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Impute missing values and compute principal components')
    parser.add_argument('--idir', default = '../dataDerived/methylBedByChrom/')    
    parser.add_argument('--odir', default = '../dataDerived/methylBedImputed/')
    parser.add_argument('--samples_file', default = '../dataDerived/masterSamplesheet.csv')
    parser.add_argument('--chrom', default = "chr18")

    parser.add_argument('--drop_MCI', action="store_true")
    parser.add_argument('--dry_run', action="store_true")
    args = parser.parse_args()

    main(args)
        
