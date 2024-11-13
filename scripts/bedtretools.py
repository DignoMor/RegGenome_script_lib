#!/usr/bin/env python

# Tools to process bedtre files outputted by PINTS

import argparse

import pandas as pd

from RGTools.BedTable import BedTable3, BedTable6

class BedTre(BedTable3):
    def __init__(self):
        self._data_df = pd.DataFrame(columns=self.column_names)
        super().__init__()
    
    @property
    def column_names(self):
        return ['chrom', 
                'start', 
                'end', 
                'name', 
                'fwdTSS', 
                'revTSS',
                ]
    
    @property
    def column_types(self):
        column_type = super().column_types
        column_type["name"] = str
        column_type["fwdTSS"] = int
        column_type["revTSS"] = int

        column_type = self._dtype2pddtype(column_type)

        return column_type

    def get_region_names(self):
        return self._data_df['name'].values
    
    def get_region_fwdTSS(self):
        return self._data_df['fwdTSS'].values
    
    def get_region_revTSS(self):
        return self._data_df['revTSS'].values
    
class BedTRETools:
    @staticmethod
    def set_parser(parser):
        subparsers = parser.add_subparsers(dest="subcommand")

        parser_get_tss_defined_boundry = subparsers.add_parser("get_tss_defined_boundry",
                                                               help="Get TSS defined boundry.",
                                                               )

        BedTRETools.set_parser_get_tss_defined_boundry(parser_get_tss_defined_boundry)

    @staticmethod
    def set_parser_get_tss_defined_boundry(parser):
        parser.add_argument("--inpath", "-I",
                            help="Input path for bedtre file (bidirectional elements).",
                            required=True,
                            dest="inpath",
                            )

        parser.add_argument("--opath", "-O",
                            help="Output path for bedtre file.",
                            required=True,
                            dest="opath",
                            )

        parser.add_argument("--tss_padding", 
                            help="Padding for TSS.",
                            required=True,
                            type=int,
                            dest="tss_padding",
                            )
        
    @staticmethod
    def main_get_tss_defined_boundry(args):
        bedtre = BedTre()
        bedtre.load_from_file(args.inpath)

        output_bt = BedTable6()
        output_df = pd.DataFrame(columns=output_bt.column_names)
        
        for region in bedtre.iter_regions():
            chrom = region["chrom"]

            if region["fwdTSS"] < region["revTSS"]: # convergent 

                start = region["revTSS"] - args.tss_padding if region["revTSS"] - args.tss_padding < region["start"] else region["start"]
                end = region["fwdTSS"] + args.tss_padding if region["fwdTSS"] + args.tss_padding + 1 > region["end"] else region["end"]
                region_type = "convergent"

            else: # divergent
                end = region["fwdTSS"] + args.tss_padding + 1
                start = region["revTSS"] - args.tss_padding
                region_type = "divergent"
                

            output_df.loc[output_df.shape[0]] = [chrom,
                                                 start, 
                                                 end, 
                                                 "_".join([region["chrom"], 
                                                           str(region["start"]),
                                                           str(region["end"]), 
                                                           region_type, 
                                                           ]), 
                                                 ".", 
                                                 ".", 
                                                 ]
        
        output_bt.load_from_dataframe(output_df)
        output_bt.write(args.opath)

    @staticmethod
    def main(args):
        if args.subcommand == "get_tss_defined_boundry":
            BedTRETools.main_get_tss_defined_boundry(args)
        else:
            raise ValueError("Invalid subcommand.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tools to process bedtre files outputted by PINTS")
    BedTRETools.set_parser(parser)
    args = parser.parse_args()
    BedTRETools.main(args)
