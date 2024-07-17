#!/usr/bin/env python

import pyBigWig
import os
import sys
import argparse
import pandas as pd

def combine_bw_intervals(rep1, rep2, chrom_length, verbose=False):
    '''
    Merge bw intervals from two replicates and return the merged intervals

    example of interval: (100, 200, 2.5)
    The interval represent a half-open interval from base 100 to 200 of signal 2.5.
    Base 100 is included in the interval but 200 is not.

    Keyword Arguments:
    rep1 - iterable of intervals from replicate 1
    rep2 - iterable of intervals from replicate 2
    chrom_length - The totoal length of the chromosome where the intervals originates
    verbose - whether to write logs to stderr
    '''

    rep1_pointer = 0
    rep2_pointer = 0

    temp_signal_stack = []
    merged_interval_list = []

    for i in range(chrom_length):
        if verbose:
            if i % 10000000 == 0:
                sys.stderr.write("processed {:d} Mb.\n".format(i//1000000))

        if rep1_pointer < len(rep1) and i >= rep1[rep1_pointer][1]:
            rep1_pointer += 1

        if rep2_pointer < len(rep2) and i >= rep2[rep2_pointer][1]:
            rep2_pointer += 1

        signal=0

        if rep1_pointer < len(rep1) and i >= rep1[rep1_pointer][0]:
            signal += rep1[rep1_pointer][2]

        if rep2_pointer < len(rep2) and i >= rep2[rep2_pointer][0]:
            signal += rep2[rep2_pointer][2]

        if signal != 0 and len(temp_signal_stack) == 0:
            temp_signal_stack.append(signal)
        elif signal != 0 and temp_signal_stack[-1] == signal:
            temp_signal_stack.append(signal)
        elif signal != 0 and temp_signal_stack[-1] != signal:
            # start new signal
            interval_start = i - len(temp_signal_stack)
            interval_end = i
            interval_sig = temp_signal_stack[-1]

            new_interval = (interval_start, interval_end, interval_sig)
            merged_interval_list.append(new_interval)

            temp_signal_stack = []
            temp_signal_stack.append(signal)
        elif len(temp_signal_stack) != 0:
            # wrap up previous signals
            interval_start = i - len(temp_signal_stack)
            interval_end = i
            interval_sig = temp_signal_stack[-1]

            new_interval = (interval_start, interval_end, interval_sig)
            merged_interval_list.append(new_interval)

            temp_signal_stack = []
        else:
            continue

    return merged_interval_list

def merge_bw_files(rep1_bw_path, rep2_bw_path, chroms, chrom_sizes,
                   opath, verbose=False):
    rep1_bw = pyBigWig.open(rep1_bw_path)
    rep2_bw = pyBigWig.open(rep2_bw_path)

    if os.path.exists(opath):
        os.remove(opath)
    merged_bw = pyBigWig.open(opath, "w")

    merged_bw.addHeader([(c, s) for c, s in zip(chroms, chrom_sizes)])

    for chrom, chrom_length in zip(chroms, chrom_sizes):

        #TODO: better fix for invalid interval error
        try:
            rep1_chrom_intervals = rep1_bw.intervals(chrom)
        except RuntimeError:
            rep1_chrom_intervals = ()
        try:
            rep2_chrom_intervals = rep2_bw.intervals(chrom)
        except RuntimeError:
            rep2_chrom_intervals = ()

        if verbose:
            sys.stderr.write("Merging {}...\n".format(chrom))

        # In case of empty chromosomes, no combining is needed
        if len(rep1_chrom_intervals) == 0 and len(rep2_chrom_intervals) == 0:
            pass
        elif len(rep1_chrom_intervals) == 0 and len(rep2_chrom_intervals) != 0:
            merged_interval_list = rep2_chrom_intervals
        elif len(rep1_chrom_intervals) != 0 and len(rep2_chrom_intervals) == 0:
            merged_interval_list = rep1_chrom_intervals
        else:
            merged_interval_list = combine_bw_intervals(rep1_chrom_intervals,
                                                        rep2_chrom_intervals,
                                                        chrom_length,
                                                        verbose=verbose)

        if len(rep1_chrom_intervals) == 0 and len(rep2_chrom_intervals) == 0:
            pass
        else:
            merged_bw.addEntries([chrom] * len(merged_interval_list),
                                [e[0] for e in merged_interval_list],
                                ends=[e[1] for e in merged_interval_list],
                                values=[e[2] for e in merged_interval_list],
                                )

    rep1_bw.close()
    rep2_bw.close()
    merged_bw.close()

def set_parser(parser):
    parser.add_argument("--rep1",
                        help="Path to the plus strand bigwig of replicate 1. (required)",
                        required=True,
                        )

    parser.add_argument("--rep2",
                        help="Path to the plus strand bigwig of replicate 2. (required)",
                        required=True,
                        )

    parser.add_argument("--opath",
                        help="Output path. (required)",
                        required=True,
                        )
    
    parser.add_argument("--chrom_size",
                        help="Path to chrom size file. (required)",
                        required=True,
                        )

    parser.add_argument("--verbose", 
                        help="Whether to write log to stderr. (False)", 
                        default=False, 
                        )


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(prog="Merge two bigwig files.")

    set_parser(arg_parser)
    args = arg_parser.parse_args()

    if args.verbose == "False":
        args.verbose = False

    chrom_size_df = pd.read_csv(args.chrom_size,
                                sep="\t",
                                names=['chr', 'size'],
                                )

    merge_bw_files(args.rep1,
                   args.rep2,
                   chroms=chrom_size_df['chr'].values,
                   chrom_sizes=chrom_size_df['size'].values,
                   opath=args.opath, 
                   verbose=args.verbose,
                   )

