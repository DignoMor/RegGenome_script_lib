#!/usr/bin/env python

import os
import re
import sys
import pyBigWig
import argparse
import numpy as np
import pandas as pd

from RGTools.utils import str2bool
from RGTools.BedTable import BedTable3, BedTable6, BedTable6Plus
from RGTools.GenomicElements import GenomicElements
from RGTools.BwTrack import BwTrack

class CountBwSig:
    @staticmethod
    def get_region_id_types():
        '''
        Return the list of region id types.
        '''
        return ["chrom_start_end", "name", "chrom_start_end_name"]

    @staticmethod 
    def get_region_id_converter(region_id_type):
        '''
        Return a function that converts a region series to a region id 
        based on the id type.

        Keyword arguments:
        - region_id_type: type of the region id
        '''
        if region_id_type == "chrom_start_end":
            return lambda x: x["chrom"] + "_" + str(x["start"]) + "_" + str(x["end"])
        elif region_id_type == "name":
            return lambda x: x["name"]
        elif region_id_type == "chrom_start_end_name":
            return lambda x: x["chrom"] + "_" + str(x["start"]) + "_" + str(x["end"]) + "_" + x["name"]
        else:
            raise Exception("Unsupported region id type ({}).".format(region_id_type))

    
    @staticmethod
    def set_parser(parser):
        GenomicElements.set_parser_genomic_element_region(parser)
        parser.add_argument("--job_name", 
                            help="Name of the job (also output file prefix).", 
                            required=True,
                            type=str, 
                            dest="job_name", 
                            )

        parser.add_argument("--sample", 
                            help="Sample name.", 
                            action="append", 
                            required=True, 
                            type=str, 
                            dest="sample_names", 
                            )

        parser.add_argument("--bw_pl", 
                            help="bigwig file for plus strand (must match sample order).", 
                            action="append", 
                            required=True, 
                            type=str, 
                            dest="bw_pls", 
                            )

        parser.add_argument("--bw_mn", 
                            help="bigwig file for minus strand (must match sample order).", 
                            action="append", 
                            type=str, 
                            dest="bw_mns", 
                            )
        
        parser.add_argument("--single_bw", 
                            help="If there is only plus strand bigwig file. If True, bw_mn will be ignored "
                                "and ignored_strandness will be set to True.", 
                            type=str2bool, 
                            dest="single_bw", 
                            default=False,
                            )

        parser.add_argument("--opath", 
                            help="Output path.", 
                            type=str, 
                            required=True, 
                            dest="opath", 
                            )
        
        parser.add_argument("--ignore_strandness", 
                            help="If to ignore strandness.", 
                            type=str2bool, 
                            default=False, 
                            dest="ignore_strandness", 
                            )
        
        parser.add_argument("--l_pad",
                            help="Padding for the left side of the region. [0]"
                                "Positive values will expand the region, "
                                "negative values will shrink the region.",
                            type=int,
                            default=0,
                            dest="l_pad",
                            )
        
        parser.add_argument("--r_pad",
                            help="Padding for the right side of the region. [0]"
                                "Positive values will expand the region, "
                                "negative values will shrink the region.",
                            type=int,
                            default=0,
                            dest="r_pad",
                            )
        
        parser.add_argument("--min_len_after_padding",
                            help="Minimum length of the region after padding. [0]", 
                            type=int, 
                            default=50, 
                            dest="min_len_after_padding", 
                            )
        
        parser.add_argument("--method_resolving_invalid_padding",
                            help="Method to resolve invalid padding. [raise]"
                                "Options: "
                                "- raise: raise an exception, "
                                "- fallback: give up on padding."
                                "- drop: drop the region.", 
                            type=str, 
                            default="raise", 
                            dest="method_resolving_invalid_padding", 
                            )
        
        parser.add_argument("--output_type",
                            help="What information is outputted. "
                                "Options: [{}].".format(", ".join(BwTrack.get_supported_quantification_type())),
                            type=str,
                            default="raw_count",
                            dest="output_type",
                            )
        
        parser.add_argument("--region_id_type", 
                            help="Type of the region id. [chrom_start_end] "
                                 "(Options: {})".format(", ".join(CountBwSig.get_region_id_types())),
                            type=str,
                            default="chrom_start_end",
                            dest="region_id_type",
                            )

    @staticmethod
    def write_output_table(count_df, region_df, output_type, opath, job_name):
        '''
        Write the output table.
        '''
        if output_type == "RPK":
            count_df.to_csv(os.path.join(opath, job_name + ".RPK.csv"))
        elif output_type == "raw_count":
            count_df.to_csv(os.path.join(opath, job_name + ".count.csv"))
        elif output_type == "PausingIndex":
            count_df.to_csv(os.path.join(opath, job_name + ".PI.csv"))

        region_df.to_csv(os.path.join(opath, job_name + ".region_info.csv"))

    @staticmethod
    def args_check_and_preprocessing(args):
        '''
        Check the validity of the input arguments and preprocess them.
        
        Returns the preprocessed arguments.
        Argument list:
        - sample_names: list of sample names
        - bw_pls: list of paths to plus strand bigwig files
        - bw_mns: list of paths to minus strand bigwig files
        - region_file_path: path to the region file
        - region_file_type: type of the region file
        - opath: output path
        - ignore_strandness: if to ignore strandness
        - region_padding: padding for each region (string format).
        - output_type: what information is outputted
        - l_pad: left padding
        - r_pad: right padding
        '''
        if args.single_bw:
            # set bw_mns the same as bw_pls for input consistency
            args.bw_mns = args.bw_pls

            # ignore strandness
            if not args.ignore_strandness:
                raise Exception("Strandness is not ignored while single_bw is set to True.")

        if len(args.sample_names) != len(args.bw_pls) or len(args.sample_names) != len(args.bw_mns):
            raise Exception("Number of samples do not match the number of bigwig files.")
        
        if not os.path.exists(args.region_file_path):
            raise Exception("Region file does not exist.")
        
        if not os.path.exists(args.opath):
            os.makedirs(args.opath)

        if not args.region_file_type in GenomicElements.get_region_file_suffix2class_dict().keys():
            raise Exception("Unsupported region file type ({}).".format(args.region_file_type))

        if not args.output_type in BwTrack.get_supported_quantification_type():
            raise Exception("Unsupported output type ({}).".format(args.output_type))

        if not args.min_len_after_padding >= 1:
            raise Exception("Minimum length after padding should be positive.")

        if not args.method_resolving_invalid_padding in ["raise", "fallback", "drop"]:
            raise Exception("Unsupported method to resolve invalid padding ({}).".format(args.method_resolving_invalid_padding))

        if args.output_type == "PausingIndex" and (args.l_pad < 0 or args.r_pad < 0):
            sys.stderr.write("Calculating pausing index with negative padding is not recommended.\n")
        
        if not args.region_id_type in CountBwSig.get_region_id_types():
            raise Exception("Unsupported region id type ({}).".format(args.region_id_type))

        return args

    def parse_region_input(region_file, file_type, region_id_type):
        '''
        Parse the input file and return a dataframe include region info.
        The returned df could have type-specific additional information that can be 
        used in further analysis, but the followings are always included:
        - "chrom"
        - "start"
        - "end"
        - "name"
        - "score"
        - "strand"
        
        The start and end loc will following the bed file convention (0-base, half-open).
        The index of the region df will be the unique identifier for the region 
        in the output matrix.
        '''
        genomic_elements = GenomicElements(region_path=region_file,
                                           genome_path=None,
                                           region_file_type=file_type,
                                           )

        region_bed_table=genomic_elements.get_region_bed_table()
        
        if file_type == "bed3":
            new_region_bed_table = BedTable6()
            new_region_bed_table.load_from_BedTable3(region_bed_table)
            region_bed_table = new_region_bed_table

        #TODO: Switch from df to BedTables
        region_df = region_bed_table.to_dataframe()
        region_df.fillna(".", inplace=True)
        region_df.index = region_df.agg(CountBwSig.get_region_id_converter(region_id_type), axis=1)

        return region_df

    @staticmethod
    def main(args):
        args = CountBwSig.args_check_and_preprocessing(args)

        region_df = CountBwSig.parse_region_input(args.region_file_path, 
                                                  args.region_file_type, 
                                                  args.region_id_type, 
                                                  )

        if args.ignore_strandness:
            region_df["strand"] = "."

        count_df = pd.DataFrame(index=region_df.index, 
                                columns=args.sample_names, 
                                dtype=int, 
                                )

        for sample ,bw_pl_path, bw_mn_path in zip(args.sample_names, args.bw_pls, args.bw_mns):
            bed_track = BwTrack(bw_pl_path=bw_pl_path,
                                bw_mn_path=bw_mn_path,
                                signle_bw=args.single_bw,
                                )

            for region_id, region_info in region_df.iterrows():
                count_df.loc[region_id, sample] = bed_track.count_single_region(region_info["chrom"],
                                                                                region_info["start"],
                                                                                region_info["end"],
                                                                                region_info["strand"],
                                                                                args.output_type,
                                                                                args.l_pad,
                                                                                args.r_pad,
                                                                                args.min_len_after_padding,
                                                                                args.method_resolving_invalid_padding,
                                                                                )
            
        count_df.dropna(inplace=True,
                        axis=0,
                        )
        
        region_df = region_df.loc[count_df.index]

        CountBwSig.write_output_table(count_df, region_df, args.output_type, args.opath, args.job_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Extracting Counts from BWs")
    CountBwSig.set_parser(parser=parser)
    args = parser.parse_args()

    CountBwSig.main(args)
