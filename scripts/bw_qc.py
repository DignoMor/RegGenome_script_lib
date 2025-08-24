#! /usr/bin/env python

import argparse
import sys

import pandas as pd

from RGTools.utils import str2bool
from RGTools.BwTrack import SingleBwTrack, PairedBwTrack

class BwQc:
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--bw_pl", 
                            type=str, 
                            required=True, 
                            help="The bigwig file to be checked.", 
                            )

        parser.add_argument("--bw_mn", 
                            type=str, 
                            required=False, 
                            default=None,
                            help="The bigwig file to be checked (single bw if not provided).", 
                            )

        parser.add_argument("--chrom_sizes", 
                            type=str, 
                            required=True, 
                            help="The chromosome sizes file.", 
                            )

        parser.add_argument("--statistic", 
                            type=str, 
                            required=True, 
                            help="Statistic to be checked.", 
                            )
        
        parser.add_argument("--opath", 
                            type=str, 
                            default="stdout",
                            help="The output path of the result.", 
                            )

    @staticmethod
    def available_statistics():
        return ["total_counts"]

    @staticmethod
    def args_check_and_preprocessing(args):
        if args.statistic not in BwQc.available_statistics():
            raise ValueError(f"Invalid statistic: {args.statistic}")

        return args

    @staticmethod
    def calculate_statistics(bw_pl, bw_mn, chrom_sizes, statistic):
        chrom_size_df = pd.read_csv(chrom_sizes, sep="\t", header=None, names=["chr", "size"])

        if statistic == "total_counts":
            if not bw_mn:
                bw_track = SingleBwTrack(bw_pl)
            else:
                bw_track = PairedBwTrack(bw_pl, bw_mn)

            output_stat = 0

            for _, row in chrom_size_df.iterrows():
                chrom = row["chr"]
                chrom_size = row["size"]

                output_stat += bw_track.count_single_region(chrom, 
                                                            0, 
                                                            chrom_size, 
                                                            output_type="raw_count", 
                                                            strand=".", 
                                                            )

            return output_stat

    @staticmethod
    def main(args):
        statistcs = BwQc.calculate_statistics(args.bw_pl, args.bw_mn, args.chrom_sizes, args.statistic)

        if args.opath == "stdout":
            sys.stdout.write("{:.2f}\n".format(statistcs))
        else:
            with open(args.opath, "w") as f:
                f.write("{:.2f}".format(statistcs))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    BwQc.set_parser(parser)
    args = parser.parse_args()
    args = BwQc.args_check_and_preprocessing(args)
    BwQc.main(args)