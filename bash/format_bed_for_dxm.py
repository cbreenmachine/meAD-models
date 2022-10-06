#Input file name (REQUIRED). Format should be bed-like: chr pos pos2 regionName methylLevel coverage
import pandas as pd
import argparse

def main(args):
    df = pd.read_csv(args.ifile, sep="\t")

    #Use one region for the whole sample
    df['regionName'] = "all"
    df['methylLevel'] = df['methylated'] / df['coverage']

    #Need columns in this specific order
    df.reindex(['chrom', 'chromStart', 'chromEnd', 'regionName', 'methylLevel', 'coverage'], axis=1)
    df[df['coverage'] > 5].to_csv(args.ofile, sep="\t", index=False, header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ifile', default = "../../dataDerived/methylBedBySample/116.bed")
    parser.add_argument('--ofile', default = "../../dataDerived/methylBedForDXM/116.dxm.bed")
    args = parser.parse_args()

    main(args)
#END