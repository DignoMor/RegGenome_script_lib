#!/usr/bin/env python3

import sys
import pyBigWig
import argparse

import pandas as pd

class FilterBw:
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--inpath", "-I", 
                            help="Path to the bigwig file to filter.", 
                            required=True, 
                            type=str, 
                            )

        parser.add_argument("--opath", "-O",
                            help="Path to the output bigwig file.",
                            required=True,
                            type=str,
                            )
        
        parser.add_argument("--chrom_size",
                            help="Chrom_size file used to filter the bigwig file.",
                            required=True,
                            type=str,
                            )

    @staticmethod
    def main(args):
        input_bw = pyBigWig.open(args.inpath)
        output_bw = pyBigWig.open(args.opath, "w")

        chrom_size_df = pd.read_csv(args.chrom_size, 
                                    sep="\t", 
                                    header=None, 
                                    names=["chrom", "size"], 
                                    )

        output_bw.addHeader([(chrom, size) for chrom, size in zip(chrom_size_df["chrom"], chrom_size_df["size"])])

        for chrom in chrom_size_df["chrom"]:
            if chrom in input_bw.chroms():
                intervals = input_bw.intervals(chrom)
                if intervals:
                    output_bw.addEntries([chrom] * len(intervals),
                                        [e[0] for e in intervals],
                                        ends=[e[1] for e in intervals],
                                        values=[e[2] for e in intervals],
                                        )
            else:
                sys.stderr.write(f"WARNING: Chrom {chrom} not found in {args.inpath}\n")

        input_bw.close()
        output_bw.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter a bigwig file by chromosomes.")

    FilterBw.set_parser(parser)

    args = parser.parse_args()

    FilterBw.main(args)


