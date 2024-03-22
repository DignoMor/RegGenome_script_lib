#!/usr/bin/env python

# downstream processing of count tables

import argparse

import numpy as np
import pandas as pd

def set_parser(parser):
    parser.add_argument("--inpath", "-I", 
                        help="Input path for count table.", 
                        required=True, 
                        dest="inpath", 
                        )
    
    parser.add_argument("--cpm_opath", 
                        help="Output path for CPM table.", 
                        default=None, 
                        dest="cpm_opath", 
                        )

def count_df2cpm_df(count_df):
    CPM_table = count_df / count_df.sum(axis=0) * 1e6
    return CPM_table

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Count Table Tool.")

    set_parser(parser)

    args = parser.parse_args()

    count_df = pd.read_csv(args.inpath, 
                           index_col=0, 
                           )

    if args.cpm_opath:
        CPM_table = count_df2cpm_df(count_df)
        CPM_table.to_csv(args.cpm_opath)

