#!/usr/bin/env python

import argparse

import pandas as pd

from RGTools.BedTable import BedTable6

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

    @staticmethod
    def main(args):
        input_bt = BedTable6()
        output_bt = BedTable6()

        input_bt.load_from_file(args.inpath)
        output_df = pd.DataFrame(columns=output_bt.column_names)

        for region in input_bt.iter_regions():
            if region["strand"] == "+":
                start = region["start"] - args.upstream_pad
                end = region["end"] + args.downstream_pad
            elif region["strand"] == "-":
                start = region["start"] - args.downstream_pad
                end = region["end"] + args.upstream_pad
            else:
                raise ValueError("Invalid strand value: {}".format(region["strand"]))
            
            output_df.loc[output_df.shape[0]] = [region["chrom"], start, end, region["name"], region["score"], region["strand"]]

        output_bt.load_from_dataframe(output_df)
        
        output_bt.write(args.opath)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Add padding to a BED file")
    PaddingBed.set_parser(parser)
    
    args = parser.parse_args()

    PaddingBed.main(args)
    
