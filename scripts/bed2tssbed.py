#!/usr/bin/env python

import argparse
import sys
import re

import pandas as pd
import numpy as np

from RGTools.BedTable import BedTable6, BedTable6Plus

class Bed2TSSBED:
    @staticmethod
    def get_output_site_types():
        return ["TSS", "center"]

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
        
        parser.add_argument("--output_site", 
                            help="The site for output. [TSS] ({})".format(", ".join(Bed2TSSBED.get_output_site_types())), 
                            default="TSS",
                            type=str, 
                            )

    @staticmethod
    def sanity_check(args):
        input_bed_table = Bed2TSSBED.read_input(args.bed_in, args.region_file_type)
        # check strand information
        if args.output_site == "TSS":
            for strand in np.unique(input_bed_table.get_region_strands()):
                if strand not in ("+", "-"):
                    raise Exception("Insufficient strand information!")
        
        if not args.output_site in Bed2TSSBED.get_output_site_types():
            raise Exception("Unrecognized output site type: {}".format(args.output_site))

    def get_site_coord(region, output_site_type):
        if output_site_type == "TSS":
            if region["strand"] == "+":
                site = region["start"]
            elif region["strand"] == "-":
                site = region["end"] - 1
        elif output_site_type == "center":
            site = int((region["start"] + region["end"]) // 2)
        
        return site

    @staticmethod
    def main(args):
        Bed2TSSBED.sanity_check(args)

        input_bed_table = Bed2TSSBED.read_input(args.bed_in, args.region_file_type)
        
        output_bt = input_bed_table._clone_empty()
        output_region_list = []

        for region in input_bed_table.iter_regions():
            output_coord = Bed2TSSBED.get_site_coord(region, args.output_site)
            
            output_region_dict = region.to_dict()
            output_region_dict["start"] = output_coord
            output_region_dict["end"] = output_coord + 1
            output_region_list.append(output_region_dict)

        output_df = pd.DataFrame(output_region_list, 
                                 columns=output_bt.column_names)

        output_bt.load_from_dataframe(output_df)
        output_bt.write(args.bed_out)


if __name__ == "__main__":
    sys.stderr.write("Warning: This script is deprecated."
                     "Use genomicelement_tool.py bed2tssbed for the same functionality.\n")

    parser = argparse.ArgumentParser(prog="Gen TSS bed from bed6 files.")
    Bed2TSSBED.set_parser(parser=parser)
    args = parser.parse_args()

    Bed2TSSBED.main(args)
