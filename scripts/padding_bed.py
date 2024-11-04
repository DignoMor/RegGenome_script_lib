#!/usr/bin/env python

import argparse

import pandas as pd

from RGTools.BedTable import BedTable6, BedTable3
from RGTools.utils import str2bool

class PaddingBed:
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--inpath", 
                            help="The bed6 file path to add padding to", 
                            type=str,
                            required=True, 
                            )
        
        parser.add_argument("--upstream_pad",
                            help="Amount to extend to the upstream of the region. "
                                 "Positive value will expand the region and negative value will shrink the region.",
                            type=int,
                            required=True,
                            )
        
        parser.add_argument("--downstream_pad",
                            help="Amount to extend to the downstream of the region. "
                                 "Positive value will expand the region and negative value will shrink the region.",
                            type=int,
                            required=True,
                            )
        
        parser.add_argument("--opath",
                            help="Output path for the padded BED file",
                            type=str,
                            required=True,
                            )
        
        parser.add_argument("--input_file_type",
                            help="File type of the input region file. "
                            "If bed3 is given, all regions will be assumed "
                            "to be on the positive strand. [bed6] "
                            "(Options:{})".format(", ".join(PaddingBed.get_supported_region_file_types())), 
                            default="bed6",
                            type=str, 
                            )
        
        parser.add_argument("--ignore_strand",
                            help="Ignore the strand information in the input file. ", 
                            default=False,
                            type=str2bool,
                            )

    @staticmethod
    def args_sanity_check(args):
        if not args.input_file_type in PaddingBed.get_supported_region_file_types():
            raise ValueError("Invalid input file type: {}".format(args.input_file_type))
        
        if args.input_file_type == "bed3":
            if not args.ignore_strand:
                raise ValueError("Strand info not available in bed3 file. ")
        
        return args

    @staticmethod
    def get_supported_region_file_types():
        '''
        Return the supported region file types.
        '''
        return ["bed6", "bed3"]
    
    @staticmethod
    def read_region(region_file_path, region_file_type):
        '''
        Read a region file and return a BedTable object.
        '''
        if region_file_type == "bed6":
            bt = BedTable6()
            bt.load_from_file(region_file_path)
        elif region_file_type == "bed3":
            bt = BedTable3()
            bt.load_from_file(region_file_path)
        
        return bt

    @staticmethod
    def padding_region(region, upstream_pad, downstream_pad, 
                       ignore_strand=False, 
                       ):
        if ignore_strand:
            start = region["start"] - upstream_pad
            end = region["end"] + downstream_pad
        elif region["strand"] == "+":
            start = region["start"] - upstream_pad
            end = region["end"] + downstream_pad
        elif region["strand"] == "-":
            start = region["start"] - downstream_pad
            end = region["end"] + upstream_pad
        else:
            raise ValueError("Invalid strand value: {}".format(region["strand"]))
        
        
        return {"chrom": region["chrom"],
                "start": start,
                "end": end, 
                }

    @staticmethod
    def main(args):
        args = PaddingBed.args_sanity_check(args)

        input_bt = PaddingBed.read_region(args.inpath, args.input_file_type)
        output_bt = input_bt._clone_empty()
        output_df = input_bt.to_dataframe()

        for i, region in enumerate(input_bt.iter_regions()):
            padded_region = PaddingBed.padding_region(region, 
                                                      args.upstream_pad, 
                                                      args.downstream_pad, 
                                                      ignore_strand=args.ignore_strand,
                                                      )
            
            output_df.loc[i, "start"] = padded_region["start"]
            output_df.loc[i, "end"] = padded_region["end"]

        output_bt.load_from_dataframe(output_df)
        
        output_bt.write(args.opath)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Add padding to a BED file")
    PaddingBed.set_parser(parser)
    
    args = parser.parse_args()

    PaddingBed.main(args)
    
