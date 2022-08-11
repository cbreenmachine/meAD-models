import pandas as pd
import argparse
import os
import numpy as np
from sklearn.decomposition import IncrementalPCA

# os.chdir('/z/Comp/kelesgroup/loadmethyseq/methylome-diff-analysis/scripts-v2/')

def read_file(ifile):
    z = pd.read_csv(ifile, sep = "\t")
    return(z)

def read_and_join_all(file_list):
    '''
    Parameters
    ----------
        file_path : str
            relative path to one chromosome
    
    Returns
    -------
    '''
    print("Reading in all files")
    df_raw = pd.concat([read_file(x) for x in file_list])

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


def run_pca(X, n_components = 70):
    '''computes principal components incrementally
    '''
    #TODO: allow for normal PCA
    ipca = IncrementalPCA(n_components = n_components, batch_size=10000)
    X_ipca = ipca.fit_transform(X)
    
    pcs = pd.DataFrame(X_ipca)
    pcs.columns = ['PC'+str(x+1) for x in range(n_components)]

    data = {'var_explained': ipca.explained_variance_ratio_, 
            'PCs' : [x+1 for x in range(n_components)]}
    var = pd.DataFrame(data)

    return pcs, var


def make_name(idir, chrom, sample_num):
    return os.path.join(idir, chrom + "." + str(sample_num) + ".bed")


def main(args):

    # Handle directory structure
    odir = args.odir
    if not os.path.exists(args.odir):
        os.makedirs(args.odir)

    idir, samples_file, chrom = args.idir, args.samples_file, args.chrom

    samples_df = pd.read_csv(samples_file)


    valid_files = []
    for x in samples_df['RA_id']:
        ff = make_name(idir, chrom, x)
        if os.path.exists(ff):
            valid_files.append(ff)
        else:
            print("DNE: " + ff)    

    # This takes awhile
    df_raw = read_and_join_all(valid_files)
    df = df_raw.pivot(index='chromStart', columns='sample', values=['methylated', 'coverage'])
    
    # Split into coverage and methylated dataframes 
    # (Probably a way to do this with multi-indexing, but this is easier to understand
    # and more R-like)
    # reset_index gives chromStart as a column again, so hold off until imputation is done
    Cov = df['coverage'].copy(deep=True).rename_axis(None, axis=1)
    M = df['methylated'].copy(deep=True).rename_axis(None, axis=1)
    
    Cov.columns = Cov.columns.astype(str)
    M.columns = M.columns.astype(str)

    print("Imputing missing values...")
    for group in list(set(samples_df['diagnostic_group'])):
        print(group)
        cols = [str(x) for x in samples_df[samples_df.diagnostic_group == group]['RA_id']]
        print(str(len(cols)) + " in group " + group)
        
        impute_col_mean(M, cols)
        impute_col_mean(Cov, cols)
      
    ix_1 = np.where(Cov.isna().sum(axis=1) > 0)[0]
    ix_2 = np.where(M.isna().sum(axis=1) > 0)[0]
    ix_3 = np.where((Cov == 0).sum(axis=1) > 0)[0]

    Cov.reset_index(inplace=True)
    M.reset_index(inplace=True)  

    ix_to_drop = list(set(np.concatenate([ix_1, ix_2, ix_3])))
    print(str(len(ix_to_drop)) + " of " + str(len(Cov)) + " positions will be dropped. ")

    Cov.drop(ix_to_drop, axis=0, inplace=True)
    M.drop(ix_to_drop, axis=0, inplace=True)

    ofile = os.path.join(args.odir, chrom + ".M.bed")
    M.to_csv(ofile, sep='\t', index=False)

    ofile = os.path.join(args.odir, chrom + ".Cov.bed")
    Cov.to_csv(ofile, sep='\t', index=False)

    # Finally, run PCA
    X = M.drop('chromStart', axis=1).div(Cov.drop('chromStart', axis=1)).transpose()
    pcs, var = run_pca(X)

    sample_id = [os.path.basename(x).replace( chrom + '.', '').replace('.bed', '') for x in valid_files]
    pcs['sample_id'] = sample_id

    ofile = os.path.join(args.odir, chrom + ".PCs.bed")
    pcs.to_csv(ofile, sep='\t', index=False)

    ofile = os.path.join(args.odir, chrom + ".PCs.varexp.bed")
    var.to_csv(ofile, sep='\t', index=False)

  
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Impute missing values and compute principal components')
    parser.add_argument('--idir', default = '../data/06-counts-by-chrom/') 
    parser.add_argument('--samples_file', default = '../data/07-counts-by-chrom-imputed-subset1/master-42-subsampled-1.csv')
    parser.add_argument('--chrom', default = "chr18")
    parser.add_argument('--odir', default = '../data/07-counts-by-chrom-imputed-subset1/')
    args = parser.parse_args()

    main(args)
        
