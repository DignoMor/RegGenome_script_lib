#!/usr/bin/env python

# downstream processing of count tables

import argparse
import sys
import os

import numpy as np
import pandas as pd

import statsmodels.api as sm

from RGTools.BedTable import BedTable6, BedTable6Plus
from RGTools.utils import str2bool, str2none

class CountTableTool:
    @staticmethod
    def main(args):
        if args.subcommand == "per_million_normalization":
            CountTableTool.per_million_normalization_main(args)
        elif args.subcommand == "cat_table":
            CountTableTool.cat_table_main(args)
        elif args.subcommand == "substitute_gene_id":
            CountTableTool.substitute_gene_id_main(args)
        elif args.subcommand == "divide_table":
            CountTableTool.divide_table_main(args)
        elif args.subcommand == "compute_tissue_tstat":
            CountTableTool.compute_tissue_tstat_main(args)
        elif args.subcommand == "count_table2bed":
            CountTableTool.count_table2bed_main(args)
        elif args.subcommand == "combine_tissue_counts":
            CountTableTool.combine_tissue_counts_main(args)
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

        parser_divide_table = subparsers.add_parser("divide_table",
                                                    help="Divide operation between count tables.",
                                                    )
        
        CountTableTool.set_parser_divide_table(parser_divide_table)

        parser_compute_tissue_tstat = subparsers.add_parser("compute_tissue_tstat",
                                                            help="Compute t-statistic for tissue specificity."
                                                                 "Based on Finucane et al. 2018 Nature Genetics.",
                                                            )
        
        CountTableTool.set_parser_compute_tissue_tstat(parser_compute_tissue_tstat)

        parser_count_table2bed = subparsers.add_parser("count_table2bed",
                                                       help="Convert tstat table to bed file.",
                                                       )
        
        CountTableTool.set_parser_count_table2bed(parser_count_table2bed)

        parser_combine_tissue_counts = subparsers.add_parser("combine_tissue_counts",
                                                             help="Combine tissue counts.",
                                                             )

        CountTableTool.set_parser_combine_tissue_counts(parser_combine_tissue_counts)

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

    def set_parser_compute_tissue_tstat(parser):
        parser.add_argument("--inpath", "-I", 
                            help="Input path for count table.", 
                            required=True, 
                            dest="inpath", 
                            )
        
        parser.add_argument("--tissue_labels", 
                            help="Tissue labels for each sample. comma seperated string. "
                                 "(e.g. tissue1,tissue2,tissue3,tissue3) "
                                 "If None, column names will be used.",
                            default=None,
                            type=str2none, 
                            dest="tissue_labels", 
                            )

        parser.add_argument("--missing",
                            help="Method for handling missing values.",
                            default="raise",
                            choices=["raise", "drop"],
                            )

        parser.add_argument("--opath", 
                            help="Output path for tstat table.", 
                            default="stdout", 
                            dest="opath", 
                            )

    def set_parser_combine_tissue_counts(parser):
        CountTableTool.set_parser_compute_tissue_tstat(parser)

    def set_parser_divide_table(parser):
        parser.add_argument("--inpath_numerator",
                            help="Input paths for count tables as numerators.",
                            required=True,
                            )

        parser.add_argument("--inpath_denominator",
                            help="Input paths for count tables as denominators.",
                            required=True,
                            )
        
        parser.add_argument("--min_numerator",
                            help="Minimum value for numerator. Filtered values will be set as NA.",
                            default=0,
                            type=float,
                            )

        parser.add_argument("--min_denominator",
                            help="Minimum value for denominator. Filtered values will be set as NA.",
                            default=0,
                            type=float,
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
    def set_parser_count_table2bed(parser_count_table2bed):
        '''
        Convert tstat table to bed file.
        '''
        parser_count_table2bed.add_argument("--inpath", "-I", 
                                            help="Input path for tstat table.", 
                                            required=True, 
                                            dest="inpath", 
                                            )
        
        parser_count_table2bed.add_argument("--opath", "-O",
                                            help="Output path for bed file.", 
                                            default="stdout", 
                                            dest="opath", 
                                            )
        
        parser_count_table2bed.add_argument("--percentage",
                                            help="Percentage of top tstat to be included in the bed file.", 
                                            default=0.2, 
                                            type=float, 
                                            dest="percentage", 
                                            )
        
        parser_count_table2bed.add_argument("--region_info",
                                            help="Path to region info csv.", 
                                            required=True, 
                                            dest="region_info", 
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
    def check_column_match(df1, df2):
        if not (df1.columns == df2.columns).all():
            raise ValueError("Column mismatch.")

        return None
    
    @staticmethod
    def compute_tstat_one_elem(X, Y, missing="raise"):
        '''
        Compute t-statistics for one element.
        Adapted from You Chen's script.

        Keyword arguments:
        - X: labels for tissue type. np.array with shape (n, 1).
          n being the number of sampels. For each sample, X has value 
          1 for the tissue of interest and -1 for other tissues. 
        - Y: expression matrix. np.array with shape (n, 1).
        - missing: method handling missing values.
        '''
        X_df = pd.DataFrame(X, columns=["X1"])
        X_df = sm.add_constant(X_df["X1"])

        # Return np.nan if Y is all nan
        if np.all(np.isnan(Y)):
            return np.nan

        # Return np.nan if X is singular after removing nans
        if not np.linalg.matrix_rank(X_df.loc[~np.isnan(Y)]) > 1:
            return np.nan

        model = sm.OLS(Y, X_df, missing=missing).fit()
        t = model.tvalues["X1"]
        return t

    @staticmethod
    def per_million_normalization_main(args):
        input_df = CountTableTool.read_input_df(args.inpath)
        output_df = input_df / input_df.sum(axis=0) * 1e6
        CountTableTool.write_output_df(output_df, args.opath)
        return None
    
    @staticmethod
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

    @staticmethod
    def divide_table_main(args):
        numerator_df = CountTableTool.read_input_df(args.inpath_numerator)
        denominator_df = CountTableTool.read_input_df(args.inpath_denominator)

        CountTableTool.check_index_match(numerator_df, denominator_df)
        CountTableTool.check_column_match(numerator_df, denominator_df)

        output_df = pd.DataFrame(np.nan,
                                 index=numerator_df.index,
                                 columns=numerator_df.columns,
                                 )
        for c in output_df.columns:
            logical_pass_filter = (denominator_df[c] >= args.min_denominator) & (numerator_df[c] >= args.min_numerator)
            output_df.loc[logical_pass_filter, c] = numerator_df.loc[logical_pass_filter, c] / denominator_df.loc[logical_pass_filter, c]
        
        CountTableTool.write_output_df(output_df, args.opath)
    
    @staticmethod
    def compute_tissue_tstat_main(args):
        input_df = CountTableTool.read_input_df(args.inpath)
        if not args.tissue_labels:
            tissue_labels = np.array(input_df.columns)
        else:
            tissue_labels = np.array(args.tissue_labels.split(","))

        unique_tissues = np.unique(tissue_labels)

        output_array = np.zeros((input_df.shape[0], len(unique_tissues)))

        for tissue_ind, tissue in enumerate(unique_tissues):
            X = tissue_labels == tissue
            X = np.array([1 if x else -1 for x in X])

            for elem_ind, row in enumerate(input_df.iterrows()):
                elem = row[0]
                elem_info = row[1]
                Y = elem_info.values

                tstat = CountTableTool.compute_tstat_one_elem(X, 
                                                              Y, 
                                                              missing=args.missing, 
                                                              )
                output_array[elem_ind, tissue_ind] = tstat

        output_df = pd.DataFrame(output_array,
                                 index=input_df.index,
                                 columns=unique_tissues,
                                 )

        CountTableTool.write_output_df(output_df, args.opath)

    @staticmethod
    def combine_tissue_counts_main(args):
        input_df = CountTableTool.read_input_df(args.inpath)
        if not args.tissue_labels:
            tissue_labels = np.array(input_df.columns)
        else:
            tissue_labels = np.array(args.tissue_labels.split(","))

        unique_tissues = np.unique(tissue_labels)

        output_array = np.zeros((input_df.shape[0], len(unique_tissues)))

        for tissue_ind, tissue in enumerate(unique_tissues):
            X = tissue_labels == tissue
            X = np.array([1 if x else -1 for x in X])

            for elem_ind, row in enumerate(input_df.iterrows()):
                elem = row[0]
                elem_info = row[1]
                elem_counts = elem_info.values

                if np.all(np.isnan(elem_counts)):
                    output_array[elem_ind, tissue_ind] = np.nan
                else:
                    num_nan = np.sum(np.isnan(elem_counts))
                    num_data = len(elem_counts) - num_nan
                    count_sum = np.sum(elem_counts[~np.isnan(elem_counts)])/num_data * len(elem_counts)
                    output_array[elem_ind, tissue_ind] = count_sum

        output_df = pd.DataFrame(output_array,
                                 index=input_df.index,
                                 columns=unique_tissues,
                                 )

        CountTableTool.write_output_df(output_df, args.opath)

    @staticmethod
    def count_table2bed_main(args):
        input_df = CountTableTool.read_input_df(args.inpath)
        region_info_df = CountTableTool.read_region_info_df(args.region_info)

        if len(region_info_df.columns) == 6:
            init_output_bt = lambda: BedTable6()
        else:
            extra_columns = [c for c in region_info_df.columns if c not in ["chrom", "start", "end", "name", "score", "strand"]]
            init_output_bt = lambda: BedTable6Plus(extra_column_names=extra_columns, 
                                                   extra_column_dtype=[str] * (len(extra_columns)),
                                                   )

        for cellline in input_df.columns:
            stats = input_df[cellline].values
            stats = stats[~ np.isnan(stats)]
            cutoff = np.percentile(stats, 100 - args.percentage * 100)
            regions_to_keep = input_df[input_df[cellline] > cutoff].index

            output_bt = init_output_bt()
            output_bt.load_from_dataframe(region_info_df.loc[regions_to_keep, :])
            output_bt.write(os.path.join(args.opath, 
                                         cellline + ".top.{:d}%.bed".format(int(args.percentage * 100)), 
                                         ))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Count Table Tool.")

    CountTableTool.set_parser(parser)
    args = parser.parse_args()
    CountTableTool.main(args)
