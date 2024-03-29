#!/usr/bin/env python

import argparse
import re

import pandas as pd

def set_parser(parser):
    parser.add_argument("--bed_in", "-I",
                        help="Path to input bed6 file.",
                        required=True, 
                        )

    parser.add_argument("--bed_out", '-O',
                        help="output path for bed results.",
                        required=True, 
                        )

    parser.add_argument("--window_size", "-w",
                        help="The window size [250-250]. The default is to return 501 bp window with tss in middle.", 
                        default="250-250")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Gen TSS bed from bed6 files.")
    set_parser(parser=parser)
    args = parser.parse_args()

    bed_df = pd.read_csv(args.bed_in, 
                         sep="\t", 
                         names=["chr_name", 
                                "start", 
                                "end", 
                                "gene_id", 
                                "score", 
                                "strand", 
                                ], 
                         )
    
    # check strand information
    for strand in bed_df["strand"].unique():
        if strand not in ("+", "-"):
            raise Exception("Insufficient strand information!")

    match = re.search("([0-9]*)-([0-9]*)", args.window_size)
    if not match:
        raise Exception("Unrecognized window size: {}".format(args.window_size))

    l_pad = int(match.group(1))
    r_pad = int(match.group(2))
    
    tss_list = bed_df["start"].values
    tss_list[bed_df["strand"]=="-"] = bed_df.loc[bed_df["strand"]=="-", "end"] - 1

    bed_df["start"] = [t - l_pad for t in tss_list] # bed is 0-based
    bed_df["end"] = [t + r_pad + 1 for t in tss_list]

    bed_df.to_csv(args.bed_out,
                  sep="\t",
                  header=False,
                  index=False,
                  )
