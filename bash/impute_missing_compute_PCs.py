import pandas as pd
import argparse
import os
import numpy as np
from sklearn.decomposition import IncrementalPCA
import logging



# os.chdir("/z/Comp/kelesgroup/loadmethyseq/methylome-diff-analysis/bash")

#TODO: pandas not recognizing these if you use df.dtypes
# col_types = {"chrom": "str", "chromStart": "int32", "strand": "str", "methylated": np.int8, "coverage": "int8", "sample": "int8"}

def construct_full_index(all_files):
    '''Find all unique positions so we don't have to pivot'''
    tmp = pd.concat((pd.read_csv(ff, sep="\t", usecols=["chromStart"]) for ff in all_files), axis=0)
    full_index = tmp['chromStart'].unique()
    return full_index


def read_one_file(ff, full_index):
    '''Pull in one file with chrom chromStart...methylated, coverage, sample'''
    df = pd.read_csv(ff, sep="\t").set_index("chromStart").reindex(index=full_index)

    # For example, get '100' from ./path/chr22.100.bed
    sample = os.path.basename(ff).split('.')[1]

    # column becomes sample_id (e.g. '100')
    M = df['methylated'].rename(sample)
    Cov = df['coverage'].rename(sample)

    return(M, Cov)

def read_and_join_all(all_files, full_index):
    '''Read, expand, and concatanate all M and Cov matrices'''

    # zip(*) gives you two lists instead of a list of duples
    M_list, Cov_list = zip(*[read_one_file(ff, full_index) for ff in all_files])

    # deep copy probably is no longer necessary, but doens't take that long 
    M = pd.concat(M_list, axis=1).copy(deep=True)
    Cov = pd.concat(Cov_list, axis=1).copy(deep=True)

    return M, Cov

def impute_with_mean(df, columns):
    '''Among 'columns' (samples)
    df : a pandas data frame (coverage or methylated reads)
    No smoothing because that's handled in DSS later on.
    '''
    # Adapted from https://stackoverflow.com/questions/33058590/pandas-dataframe-replacing-nan-with-row-average
    # the copy here seems to be neccessary? Avoids circular referencing?
    m = df[columns].mean(axis=1, skipna=True).round().copy(deep=True)
    for c in columns:
        df[c] = df[c].fillna(m)
    return df


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
    '''returns something like inDirectory/chr18.100.bed'''
    return os.path.join(idir, chrom + "." + str(sample_num) + ".bed")


def get_samples_by_group(df, group):
    all_group_samples = set([str(x) for x in df[df.diagnostic_group == group]['sample_id']])
    return list(all_group_samples)

def main(args):

    idir, samples_file, chrom = args.idir, args.samples_file, args.chrom
    samples_df = pd.read_csv(samples_file)

    samples_by_group = dict()

    for gg in ["CONTROL", "LOAD"]:
        samples_by_group[gg] = get_samples_by_group(samples_df, gg)

    # Handle the MCI versus no MCI logic
    if args.include_MCI:
        odir = args.odir + "controlMCILOAD"
        samples_by_group["MCI"] = get_samples_by_group(samples_df, "MCI")
    else:
        odir = args.odir + "controlLOAD"
        
    # Handle directory structure
    if not os.path.exists(odir):
        os.makedirs(odir)

    all_files = []
    l = list(samples_by_group.values())

    # This combines with the dictionary logic should only allow us
    # to load in CONTROL/LOAD if we don't specify we want MCI
    # https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
    samples_of_interest = [item for sublist in l for item in sublist]

    for x in samples_of_interest:
        ff = make_name(idir, chrom, x)
        if os.path.exists(ff):
            all_files.append(ff)
    
    # reset_index gives chromStart as a column again, so hold off until imputation is done
    full_index = construct_full_index(all_files)
    print("Constructed full index")

    M, Cov = read_and_join_all(all_files, full_index)
    loaded_samples = set(M.columns.to_list())
    
    # DO NOT CHANGE; the underlying optimization problem 
    # breaks down when you approach 100% variance explained 
    # I.E. number of components == number of samples breaks
    n_components = int(len(loaded_samples) / 4)

    print("Imputing missing values...")
    for group in samples_by_group:
        valid_cols = list(loaded_samples.intersection(samples_by_group[group]))
        print("Group {0} has {1} samples".format(group, len(valid_cols)))

        M = impute_with_mean(M, valid_cols)
        Cov = impute_with_mean(Cov, valid_cols)
      
    print("Determining drop positions")
    ix_1 = np.where(Cov.isna().sum(axis=1) > 0)[0]
    ix_2 = np.where(M.isna().sum(axis=1) > 0)[0]
    ix_3 = np.where((Cov == 0).sum(axis=1) > 0)[0]

    # Turns chromStart into columns and index becomes 0,1,2,...
    print("Resetting indices")
    Cov.reset_index(inplace=True)
    M.reset_index(inplace=True)  

    ix_to_drop = list(set(np.concatenate([ix_1, ix_2, ix_3])))
    print(str(len(ix_to_drop)) + " of " + str(len(Cov)) + " positions will be dropped. ")

    print("Dropping 'bad' positions")
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
    parser.add_argument('--odir_root', default = '../dataDerived/methylBedImputed-')
    parser.add_argument('--samples_file', default = '../dataDerived/masterSamplesheet.csv')
    parser.add_argument('--chrom', default = "chr20")
    parser.add_argument('--include_MCI', action="store_true")
    args = parser.parse_args()

    main(args)
        
