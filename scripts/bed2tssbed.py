#!/usr/bin/env python

import argparse
import re

import pandas as pd
import numpy as np

from RGTools.BedTable import BedTable6

class Bed2TSSBED:

    @staticmethod
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

    @staticmethod
    def main(args):
        input_bed_table = BedTable6()
        input_bed_table.load_from_file(args.bed_in)
        
        # check strand information
        for strand in np.unique(input_bed_table.get_region_strands()):
            if strand not in ("+", "-"):
                raise Exception("Insufficient strand information!")

        match = re.search("([0-9]*)-([0-9]*)", args.window_size)
        if not match:
            raise Exception("Unrecognized window size: {}".format(args.window_size))

        l_pad = int(match.group(1))
        r_pad = int(match.group(2))
        
        tss_list = input_bed_table.get_start_locs()
        tss_list[input_bed_table.get_region_strands()=="-"] = input_bed_table.apply_logical_filter(input_bed_table.get_region_strands()=="-").get_end_locs() - 1

        output_bed_df = pd.DataFrame(columns=["chrom", "start", "end", "name", "score", "strand"])
        output_bed_df["chrom"] = input_bed_table.get_chrom_names()
        output_bed_df["start"] = [t - l_pad for t in tss_list] # bed is 0-based
        output_bed_df["end"] = [t + r_pad + 1 for t in tss_list]
        output_bed_df["name"] = input_bed_table.get_region_names()
        output_bed_df["score"] = input_bed_table.get_region_scores()  
        output_bed_df["strand"] = input_bed_table.get_region_strands()

        output_bed_table = BedTable6()
        output_bed_table.load_from_dataframe(output_bed_df)
        output_bed_table.write(args.bed_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Gen TSS bed from bed6 files.")
    Bed2TSSBED.set_parser(parser=parser)
    args = parser.parse_args()

    Bed2TSSBED.main(args)
