import numpy as np
import pandas as pd
import argparse


# Similar to chrom 22
N_rows = 500000
N_cols = 100

d = pd.DataFrame(np.zeros((N_rows, N_cols)))
d.to_csv("test.tsv", sep="\t")



