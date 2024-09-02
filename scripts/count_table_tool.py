#!/usr/bin/env python

# downstream processing of count tables

import argparse
import sys

import numpy as np
import pandas as pd

from RGTools.utils import str2bool

class CountTableTool:
    @staticmethod
    def main(args):
        if args.subcommand == "per_million_normalization":
            CountTableTool.per_million_normalization_main(args)
        elif args.subcommand == "cat_table":
            CountTableTool.cat_table_main(args)
        elif args.subcommand == "substitute_gene_id":
            CountTableTool.substitute_gene_id_main(args)
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

        parser_substitute_gene_id = subparsers.add_parser("substitute_gene_id",
                                                          help="Substitute gene ID.",
                                                          )
        
        CountTableTool.set_parser_substitute_gene_id(parser_substitute_gene_id)

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
    def set_parser_substitute_gene_id(parser):
        parser.add_argument("--inpath", "-I", 
                            help="Input path for count table.", 
                            required=True, 
                            dest="inpath", 
                            )
        
        parser.add_argument("--opath", "-O", 
                            help="Output path.", 
                            default="stdout", 
                            dest="opath", 
                            )

        parser.add_argument("--region_info_path", 
                            help="Path to region info csv.", 
                            required=True, 
                            dest="region_info_path", 
                            )
        
        parser.add_argument("--gene_id_col",
                            help="Column name for gene ID.", 
                            required=True, 
                            dest="gene_id_col", 
                            )
        
        parser.add_argument("--sort", 
                            help="Sort by new index.", 
                            dest="sort",
                            default=False,
                            type=str2bool,
                            )

    @staticmethod
    def read_input_df(input_path):
        return pd.read_csv(input_path, 
                           index_col=0, 
                           )

    @staticmethod
    def read_region_info_df(region_info_path):
        return pd.read_csv(region_info_path, 
                           index_col=0, 
                           )

    @staticmethod
    def write_output_df(output_df, opath):
        if opath == "stdout":
            output_df.to_csv(sys.stdout)
        else:
            output_df.to_csv(opath)

    @staticmethod
    def check_index_match(df1, df2):
        if not (df1.index == df2.index).all():
            raise ValueError("Index mismatch.")

        return None

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
    
    @staticmethod
    def substitute_gene_id_main(args):
        input_df = CountTableTool.read_input_df(args.inpath)
        region_info_df = CountTableTool.read_region_info_df(args.region_info_path)

        CountTableTool.check_index_match(input_df, region_info_df)

        new_ids = region_info_df[args.gene_id_col].values
        output_df = input_df.set_index(new_ids, 
                                       drop=True, 
                                       inplace=False, 
                                       )
        
        if args.sort:
            output_df = output_df.sort_index()

        CountTableTool.write_output_df(output_df, args.opath)

        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Count Table Tool.")

    CountTableTool.set_parser(parser)
    args = parser.parse_args()
    CountTableTool.main(args)
