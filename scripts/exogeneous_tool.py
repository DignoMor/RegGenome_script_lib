#!/usr/bin/env python3

import argparse
import sys
import re

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from Bio import SeqIO

from RGTools.BedTable import BedTable3

class ExogeneousTool:
    @staticmethod
    def set_parser(parser):
        subparsers = parser.add_subparsers(dest="subcommand")

        parser_metaplot = subparsers.add_parser("metaplot",
                                                help="Generate metaplot out of exogeneous sequence "
                                                     "with signal tracks.",
                                                )
        ExogeneousTool.set_parser_metaplot(parser_metaplot)

        parser_filter = subparsers.add_parser("filter",
                                              help="Filter exogeneous sequence by a regular expression.",
                                              )
        ExogeneousTool.set_parser_filter(parser_filter)

    @staticmethod
    def set_parser_metaplot(parser):
        parser.add_argument("--inpath", "-I",
                            help="Input path for exogeneous sequence (fasta).",
                            required=True,
                            )
        
        parser.add_argument("--region_path", 
                            help="Path to the region file for plotting. (bed3)",
                            required=True,
                            )

        parser.add_argument("--outpath", "-O",
                            help="Output path for metaplot.",
                            required=True,
                            )

        parser.add_argument("--signal_tracks_pl", 
                            help="Plus strand signal tracks for metaplot.",
                            required=True,
                            )

        parser.add_argument("--signal_tracks_mn", 
                            help="Minus strand signal tracks for metaplot.",
                            required=True,
                            )
        
        parser.add_argument("--title",
                            help="Title for the metaplot.",
                            default="Metaplot",
                            type=str, 
                            )

    @staticmethod
    def set_parser_filter(parser):
        parser.add_argument("--inpath", "-I",
                            help="Input path for exogeneous sequence (fasta).",
                            required=True,
                            )

        parser.add_argument("--outpath", "-O",
                            help="Output path for filtered exogeneous sequence.",
                            required=True,
                            )

        parser.add_argument("--regex", 
                            help="Regular expression for filtering.",
                            required=True,
                            )

    @staticmethod
    def load_fasta(fasta_path):
        '''
        Load a fasta file and return sequence ids and seqs
        
        Keyword arguments:
        - fasta_path: path to the fasta file

        Returns:
        - seq_ids: list of sequence ids
        - seqs: list of sequences
        '''
        seq_ids = []
        seqs = []

        with open(fasta_path, "r") as fasta_f:
            for record in SeqIO.parse(fasta_f, "fasta"):
                seq_ids.append(record.id)
                seqs.append(str(record.seq))

        return seq_ids, seqs
    
    @staticmethod
    def check_signal_track_arr(signal_tracks_arr: np.array, elem_size: int, 
                               num_tracks: int):
        '''
        Check if the signal track array is valid
        
        Keyword arguments:
        - signal_tracks_arr: signal track array
        - elem_size: size of each element in the signal track array
        - num_tracks: number of tracks in the signal track array

        Returns:
        - None
        '''
        if not isinstance(signal_tracks_arr, np.ndarray):
            raise ValueError("Signal track array should be numpy array.")

        if signal_tracks_arr.ndim != 2:
            raise ValueError("Signal track array should be 2D numpy array.")

        if signal_tracks_arr.shape[1] != elem_size:
            raise ValueError("Signal track array should have the same size of elements.")
        
        if signal_tracks_arr.shape[0] != num_tracks:
            raise ValueError("Signal track array should have the same number of tracks.")
        
    @staticmethod
    def plot_track(track_arr, ax):
        '''
        Plot a track on the region
        
        Keyword arguments:
        - track_arr: track array

        Returns:
        - None
        '''
        track2plot = track_arr.mean(axis=0)
        elem_size = len(track2plot)

        for i in range(0, elem_size):
            plot_ind = i - elem_size// 2
            if track2plot[i] == 0:
                continue

            if track2plot[i] < 0:
                rect = patches.Rectangle(
                    (plot_ind - 0.5, track2plot[i]), 
                    1, 
                    -track2plot[i], 
                    color="k", 
                )
            else:
                rect = patches.Rectangle(
                    (plot_ind  - 0.5, 0), 
                    1, 
                    track2plot[i], 
                    color="k", 
                )

            ax.add_patch(rect)

        
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks([])

        ax.set_xlim(-elem_size//2, elem_size//2)
        y_min = track2plot.min() if track2plot.min() < 0 else 0
        y_max = track2plot.max() if track2plot.max() > 0 else 0

        ax.set_ylim(y_min, y_max)

    @staticmethod
    def metaplot_main(args):
        region_bt = BedTable3()
        region_bt.load_from_file(args.region_path)

        region_sizes = region_bt.get_end_locs() - region_bt.get_start_locs()
        elem_size = region_sizes[0]
        num_elem = len(region_bt)

        if not (region_sizes == elem_size).all():
            raise ValueError("Elem sizes should be the same.")

        pl_tracks_arr = np.load(args.signal_tracks_pl)
        mn_tracks_arr = np.load(args.signal_tracks_mn)

        ExogeneousTool.check_signal_track_arr(pl_tracks_arr, elem_size, num_elem)
        ExogeneousTool.check_signal_track_arr(mn_tracks_arr, elem_size, num_elem)

        fig, axes = plt.subplots(2, 1)

        ExogeneousTool.plot_track(pl_tracks_arr, axes[0])
        ExogeneousTool.plot_track(mn_tracks_arr, axes[1])

        fig.suptitle(args.title)
        axes[1].set_xticks(range(-elem_size//2, elem_size//2+1, 100))
        fig.tight_layout()
        fig.savefig(args.outpath)

    @staticmethod
    def filter_main(args):
        with open(args.inpath, "r") as input_f:
            for record in SeqIO.parse(input_f, "fasta"):
                if re.match(args.regex, record.id):
                    with open(args.outpath, "a") as output_f:
                        SeqIO.write(record, output_f, "fasta")

    @staticmethod
    def main(args):
        if args.subcommand == "metaplot":
            ExogeneousTool.metaplot_main(args)

        elif args.subcommand == "filter":
            ExogeneousTool.filter_main(args)

        else:
            raise NotImplementedError


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Exogeneous Tool')
    ExogeneousTool.set_parser(parser)
    args = parser.parse_args()
    sys.exit(ExogeneousTool.main(args))

