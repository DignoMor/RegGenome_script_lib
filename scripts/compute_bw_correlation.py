#!/usr/bin/env python

# Author: Xiuqi Pan (xp76)
# date: 6/7/2023
# email: xp76@cornell.edu

import pyBigWig
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import pearsonr

def get_windowed_sums(bw, window_size, chrom, chrom_size):
    '''
    Get a numpy array of mean signals from a bw
    with fixed window size.
    '''
    sum_signals = []
    for bin_index in range(chrom_size//window_size - 1):
        sum_signal = bw.stats(chrom,
                              window_size*(bin_index),
                              window_size*(bin_index+1),
                              type="sum",
                              )
        if sum_signal[0]:
            sum_signals.append(sum_signal[0])
        else:
            sum_signals.append(0)

    return np.array(sum_signals)

def get_windowed_correlation(bw1_pl_path, bw2_pl_path, bw1_mn_path=None,
                             bw2_mn_path=None, window_size=50,
                             chroms=["chr1", "chr2"], chrom_sizes=["1000", "1000"],
                             ax=None, verbose=False):
    '''
    Get windowed sum of signals and compute correlation

    Keyword Arguments:
    bw1_pl_path - path to the first plus strand bigwig file
    bw1_mn_path - path to the first minus strand bigwig file
    bw2_pl_path - path to the second plus strand bigwig file
    bw2_mn_path - path to the second minus strand bigwig file
    window_size - window size to work with
    ax - ax to plot a scatter plot
    '''
    bw1_pl = pyBigWig.open(bw1_pl_path)
    bw2_pl = pyBigWig.open(bw2_pl_path)

    mn_used = bw1_mn_path and bw2_mn_path
    if mn_used:
        bw1_mn = pyBigWig.open(bw1_mn_path)
        bw2_mn = pyBigWig.open(bw2_mn_path)

    summed_signal1_array_list = []
    summed_signal2_array_list = []

    # loop through pls
    if verbose:
        sys.stderr.write("Processing Plus Strand signal......\n")
    for chrom, chrom_size in zip(chroms, chrom_sizes):
        if verbose:
            sys.stderr.write("Now on chromosome: {}\n".format(chrom))
        summed_signal1_chrom = get_windowed_sums(bw1_pl, window_size, chrom, chrom_size)
        summed_signal2_chrom = get_windowed_sums(bw2_pl, window_size, chrom, chrom_size)

        nonzero_filter = (summed_signal1_chrom > 0) & (summed_signal2_chrom > 0)

        summed_signal1_array_list.append(summed_signal1_chrom[nonzero_filter])
        summed_signal2_array_list.append(summed_signal2_chrom[nonzero_filter])

    # loop through mns
    if mn_used:
        for chrom, chrom_size in zip(chroms, chrom_sizes):
            summed_signal1_chrom = get_windowed_sums(bw1_mn, window_size, chrom, chrom_size)
            summed_signal2_chrom = get_windowed_sums(bw2_mn, window_size, chrom, chrom_size)

            nonzero_filter = (summed_signal1_chrom < 0) & (summed_signal2_chrom < 0)

            summed_signal1_array_list.append(-summed_signal1_chrom[nonzero_filter])
            summed_signal2_array_list.append(-summed_signal2_chrom[nonzero_filter])

    summed_signal1 = np.concatenate(summed_signal1_array_list)
    summed_signal2 = np.concatenate(summed_signal2_array_list)

    corr, _ = pearsonr(np.log10(summed_signal1), np.log10(summed_signal2))

    if ax:
        ax.scatter(np.log10(summed_signal1), np.log10(summed_signal2), s=0.1)
        ax.set_xlabel("rep1")
        ax.set_ylabel("rep2")

    bw1_pl.close()
    bw2_pl.close()
    if mn_used:
        bw1_mn.close()
        bw2_mn.close()

    return corr

def set_parser(parser):
    parser.add_argument("--job_name",
                        help="Name for the job. (required)",
                        required=True,
                        )

    parser.add_argument("--rep1_pl",
                        help="Path to the plus strand bigwig of replicate 1. (required)",
                        required=True,
                        )

    parser.add_argument("--rep2_pl",
                        help="Path to the plus strand bigwig of replicate 2. (required)",
                        required=True,
                        )

    parser.add_argument("--chrom_size",
                        help="Path to chrom size file. (required)",
                        required=True,
                        )

    parser.add_argument("--rep1_mn",
                        help="Path to the minus strand bigwig of replicate 1. [None]",
                        default=None,
                        )

    parser.add_argument("--rep2_mn",
                        help="Path to the minus strand bigwig of replicate 2. [None]",
                        default=None,
                        )

    parser.add_argument("--window_size",
                        help="Size of the window used to compute size. [50]",
                        default=50,
                        )

    parser.add_argument("--fig_opath",
                        help="Path to the output scatter figure. None if not outputed. [None]",
                        default=None,
                        )


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(prog="compute bw correlation")

    set_parser(arg_parser)
    args = arg_parser.parse_args()

    chrom_size_df = pd.read_csv(args.chrom_size,
                                sep="\t",
                                names=["chr", "size"],
                                )

    if args.fig_opath:
        fig, ax = plt.subplots(1, 1,
                               figsize=(8,8),
                               )
    else:
        ax=None

    corr = get_windowed_correlation(args.rep1_pl,
                                    args.rep2_pl,
                                    bw1_mn_path=args.rep1_mn,
                                    bw2_mn_path=args.rep2_mn,
                                    window_size=int(args.window_size),
                                    chroms=chrom_size_df["chr"].values,
                                    chrom_sizes=chrom_size_df["size"].values,
                                    ax=ax,
                                    )

    if ax:
        fig.suptitle(args.job_name+ " (corr={:.4f})".format(corr))

        if os.path.exists(args.fig_opath):
            os.remove(args.fig_opath)

        fig.savefig(args.fig_opath)

    sys.stdout.write("{:.4f}".format(corr))
