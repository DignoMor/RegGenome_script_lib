#!/usr/bin/env python

import argparse
import pyBigWig
import sys

import numpy as np
import pandas as pd


class SampleBw:
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--inpath", "-I", 
                            help="Path to the bigwig file to sample.", 
                            required=True, 
                            type=str, 
                            )

        parser.add_argument("--sample_rate", "-S", 
                            help="Sample rate to sample the bigwig file.", 
                            required=True, 
                            type=float, 
                            )

        parser.add_argument("--chrom_size", "-C", 
                            help="Path to the chromosome size file.", 
                            required=True, 
                            type=str, 
                            )

        parser.add_argument("--seed", 
                            help="Seed for the random number generator.", 
                            required=True, 
                            type=int, 
                            )

        parser.add_argument("--opath", "-O", 
                            help="Path to the output sampled bigwig file.", 
                            required=True, 
                            type=str, 
                            )

    @staticmethod
    def main(args):
        np.random.seed(args.seed)
        input_bw = pyBigWig.open(args.inpath)
        output_bw = pyBigWig.open(args.opath, "w")

        chrom_size_df = pd.read_csv(args.chrom_size, 
                                    sep="\t", 
                                    header=None, 
                                    names=["chrom", "size"], 
                                    )

        output_bw.addHeader([(chrom, size) for chrom, size in zip(chrom_size_df["chrom"], chrom_size_df["size"])])

        total_counts = 0

        for chrom in input_bw.chroms():
            total_counts += np.sum([e[2] for e in input_bw.intervals(chrom)])

        for chrom in chrom_size_df["chrom"]:
            if chrom in input_bw.chroms():
                intervals = input_bw.intervals(chrom)
                starts = []
                ends = []
                values = []
                for interval in intervals:
                    count_frac = interval[2] / total_counts
                    new_count = np.random.binomial(total_counts * args.sample_rate, count_frac)

                    if not (new_count == 0):
                        starts.append(interval[0])
                        ends.append(interval[1])
                        values.append(float(new_count))

                output_bw.addEntries([chrom] * len(starts), 
                                      starts, 
                                      ends=ends, 
                                      values=values,
                                      )
            else:
                sys.stderr.write(f"WARNING: Chrom {chrom} not found in {args.inpath}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    SampleBw.set_parser(parser)
    args = parser.parse_args()

    SampleBw.main(args)