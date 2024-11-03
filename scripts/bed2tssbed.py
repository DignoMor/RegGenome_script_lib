#!/usr/bin/env python

import argparse
import re

import pandas as pd
import numpy as np

from RGTools.BedTable import BedTable6, BedTable6Plus

class Bed2TSSBED:
    @staticmethod
    def BedTable6Gene():
        bt = BedTable6Plus(extra_column_names=["gene_symbol"], 
                           extra_column_dtype=[str], 
                           )
        return bt

    @staticmethod
    def get_region_file_type2class_dict():
        return {"bed6": BedTable6, 
                "bed6gene": Bed2TSSBED.BedTable6Gene}
    
    @staticmethod
    def read_input(input_file, region_file_type):
        if region_file_type not in Bed2TSSBED.get_region_file_type2class_dict().keys():
            raise Exception("Unrecognized region file type: {}".format(region_file_type))
        
        region_table = Bed2TSSBED.get_region_file_type2class_dict()[region_file_type]()
        region_table.load_from_file(input_file)
        
        return region_table

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

        parser.add_argument("--region_file_type", 
                            help="The type of region file. [bed6] ({})".format(", ".join(Bed2TSSBED.get_region_file_type2class_dict().keys())), 
                            default="bed6",
                            type=str
                            )

    @staticmethod
    def main(args):
        input_bed_table = Bed2TSSBED.read_input(args.bed_in, args.region_file_type)
        
        # check strand information
        for strand in np.unique(input_bed_table.get_region_strands()):
            if strand not in ("+", "-"):
                raise Exception("Insufficient strand information!")
        
        output_bt = input_bed_table._clone_empty()
        output_df = pd.DataFrame(columns=output_bt.column_names)

        for region in input_bed_table.iter_regions():
            if region["strand"] == "+":
                tss = region["start"]
            elif region["strand"] == "-":
                tss = region["end"] - 1
            
            output_region_dict = region.to_dict()
            output_region_dict["start"] = tss 
            output_region_dict["end"] = tss + 1
            output_df.loc[len(output_df)] = output_region_dict

        output_bt.load_from_dataframe(output_df)
        output_bt.write(args.bed_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Gen TSS bed from bed6 files.")
    Bed2TSSBED.set_parser(parser=parser)
    args = parser.parse_args()

    Bed2TSSBED.main(args)
