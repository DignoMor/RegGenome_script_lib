#!/usr/bin/env python

# downstream processing of count tables

import argparse
import sys

import numpy as np
import pandas as pd

class CountTableTool:
    @staticmethod
    def main(args):
        if args.subcommand == "per_million_normalization":
            CountTableTool.per_million_normalization_main(args)
        elif args.subcommand == "cat_table":
            CountTableTool.cat_table_main(args)
        else:
            raise ValueError("Invalid subcommand.")

    @staticmethod
    def set_parser(parser):
        subparsers = parser.add_subparsers(dest="subcommand")

        parser_per_million_normalization = subparsers.add_parser("per_million_normalization", 
                                                                 help="Per million normalization.", 
                                                                 )

        CountTableTool.set_parser_per_million_normalization(parser_per_million_normalization)

        parser_cat_table = subparsers.add_parser("cat_table",
                                                 help="Concatenate count tables.",
                                                 )
        
        CountTableTool.set_parser_cat_table(parser_cat_table)

    @staticmethod
    def set_parser_per_million_normalization(parser):
        parser.add_argument("--inpath", "-I", 
                            help="Input path for count table.", 
                            required=True, 
                            dest="inpath", 
                            )
        
        parser.add_argument("--opath", 
                            help="Output path.", 
                            default="stdout", 
                            dest="opath", 
                            )

    @staticmethod
    def set_parser_cat_table(parser):
        parser.add_argument("--inpath", "-I", 
                            help="Input paths for count tables.", 
                            required=True, 
                            dest="inpaths", 
                            action="append",
                            )
        
        parser.add_argument("--opath", 
                            help="Output path.", 
                            default="stdout", 
                            dest="opath",
                            )

    @staticmethod
    def read_input_df(input_path):
        return pd.read_csv(input_path, 
                           index_col=0, 
                           )

    @staticmethod
    def write_output_df(output_df, opath):
        if opath == "stdout":
            output_df.to_csv(sys.stdout)
        else:
            output_df.to_csv(opath)

    @staticmethod
    def per_million_normalization_main(args):
        input_df = CountTableTool.read_input_df(args.inpath)
        output_df = input_df / input_df.sum(axis=0) * 1e6
        CountTableTool.write_output_df(output_df, args.opath)
        return None
    
    def cat_table_main(args):
        input_dfs = [CountTableTool.read_input_df(inpath) for inpath in args.inpaths]
        output_df = pd.concat(input_dfs, axis=1)
        CountTableTool.write_output_df(output_df, 
                                       args.opath, 
                                       )
        
        if not (output_df.shape[0] == input_dfs[0].shape[0]):
            raise ValueError("Index mismatch.")

        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Count Table Tool.")

    CountTableTool.set_parser(parser)
    args = parser.parse_args()
    CountTableTool.main(args)
