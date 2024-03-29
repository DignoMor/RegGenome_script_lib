# plot a heatmap based on a count table

import argparse
import seaborn

import numpy as np
import pandas as pd

from RGTools.utils import str2bool

def set_parser(parser):
    parser.add_argument("--count_table", "-I", 
                        help="path to the count table", 
                        dest="count_table", 
                        required=True, 
                        )

    parser.add_argument("--plot_path", "-O", 
                        help="Path to the output plot.", 
                        dest="plot_path", 
                        required=True, 
                        )

    parser.add_argument("--log_transform", 
                        help="If log transform the data.", 
                        dest="log_transform", 
                        default=True, 
                        type=str2bool, 
                        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Plot count table heatmap.")

    set_parser(parser)

    args = parser.parse_args()

    count_table = pd.read_csv(args.count_table, 
                              index_col=0, 
                              )
    
    if args.log_transform:
        count_table = np.log10(count_table+1)
    fig = seaborn.clustermap(count_table)

    fig.savefig(args.plot_path)

    