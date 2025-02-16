#!/usr/bin/env python3

import argparse
import json
import sys
import re
import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from Bio import SeqIO
from scipy.spatial.distance import jensenshannon

from RGTools.utils import NumpyEncoder
from RGTools.BedTable import BedTable3, BedTable6Plus
from RGTools.ExogeneousSequences import ExogeneousSequences

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

        parser_compare_mutagenesis = subparsers.add_parser("compare_mutagenesis",
                                                           help="Compare Mutagenesis Results.",
                                                           )
        
        ExogeneousTool.set_parser_compare_mutagenesis(parser_compare_mutagenesis)

        parser_compute_track_correlation = subparsers.add_parser("compute_track_correlation",
                                                                 help="Compute correlation between tracks.",
                                                                 )
        
        ExogeneousTool.set_parser_compute_track_correlation(parser_compute_track_correlation)

    @staticmethod
    def set_parser_metaplot(parser):
        ExogeneousSequences.set_parser_exogeneous_sequences(parser)

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
        
        parser.add_argument("--negate_mn_strand", 
                            help="If to negate the minus strand signal.", 
                            default=False, 
                            )
        
        parser.add_argument("--title",
                            help="Title for the metaplot.",
                            default="Metaplot",
                            type=str, 
                            )

    @staticmethod
    def set_parser_filter(parser):
        ExogeneousSequences.set_parser_exogeneous_sequences(parser)
        parser.add_argument("--fasta_out", 
                            help="Output path for filtered fasta file.",
                            required=True,
                            )
        
        parser.add_argument("--region_out", 
                            help="Output path for filtered region file.",
                            required=True,
                            )

        parser.add_argument("--regex", 
                            help="Regular expression for filtering.",
                            required=True,
                            )

    @staticmethod
    def set_parser_compare_mutagenesis(parser):
        ExogeneousSequences.set_parser_exogeneous_sequences(parser)
        parser.add_argument("--regex_ref_seqs", 
                            help="Regular expression for reference sequences. "
                                 "Must contain one capture group that captures sequence ID.",
                            required=True,
                            )
        
        parser.add_argument("--regex_mut_seqs",
                            help="Regular expression for mutated sequences. "
                                 "Must contain one capture group that captures sequence ID.",
                            required=True,
                            )
        
        parser.add_argument("--pl_track_npy", 
                            help="Numpy file for plus strand signal.", 
                            required=True,
                            )

        parser.add_argument("--mn_track_npy", 
                            help="Numpy file for minus strand signal.", 
                            required=True,
                            )
        
        parser.add_argument("--total_count_plot_path", 
                            help="Path to the total count plot.",
                            default=None, 
                            )
        
        parser.add_argument("--opath", 
                            help="Output path for the results written in json format.",
                            required=True,
                            )
        
    @staticmethod
    def set_parser_compute_track_correlation(parser):
        parser.add_argument("--pl_track1_npy", 
                            help="Numpy file for the first plus strand signal.", 
                            required=True,
                            )

        parser.add_argument("--pl_track2_npy", 
                            help="Numpy file for the second plus strand signal.", 
                            required=True,
                            )
        
        parser.add_argument("--mn_track1_npy",
                            help="Numpy file for the first minus strand signal.",
                            default=None,
                            )

        parser.add_argument("--mn_track2_npy",
                            help="Numpy file for the second minus strand signal.",
                            default=None,
                            )
        
        parser.add_argument("--jensenshannon",
                            help="Output path for the jensenshannon distance.",
                            default=None,
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
    def match_ref_mut_regex(input_bt, ref_regex, mut_regex):
        '''
        Match reference and mutated sequences with regex
        
        Keyword arguments:
        - input_bt: input bed table
        - ref_regex: regex for reference sequences
        - mut_regex: regex for mutated sequences

        Returns:
        - result_bt: result *unsorted* BedTable6Plus object, 
                     with extra columns seq_id and seq_type.
                     seq_id is the identifier that groups ref and 
                     mut seqs together, while seq_type is either
                     ref or mut.
        '''
        region_bt_w_anno = BedTable6Plus(extra_column_names=["seq_id", "seq_type"], 
                                         extra_column_dtype=[str, str],
                                         enable_sort=False, 
                                         )
        
        region_chr_list = []
        region_start_list = []
        region_end_list = []
        region_seq_id_list = []
        region_seq_type_list = []

        for region in input_bt.iter_regions():
            seq_name = region["chrom"]

            ref_match = re.match(ref_regex, seq_name)
            mut_match = re.match(mut_regex, seq_name)

            if ref_match:
                region_chr_list.append(region["chrom"])
                region_start_list.append(region["start"])
                region_end_list.append(region["end"])
                region_seq_id_list.append(ref_match.group(1))
                region_seq_type_list.append("ref")

            elif mut_match:
                region_chr_list.append(region["chrom"])
                region_start_list.append(region["start"])
                region_end_list.append(region["end"])
                region_seq_id_list.append(mut_match.group(1))
                region_seq_type_list.append("mut")

            else:
                raise ValueError("Sequence name {} does not match any regex.".format(seq_name))
        
        region_df_w_anno = pd.DataFrame({"chrom": region_chr_list,
                                         "start": region_start_list,
                                         "end": region_end_list,
                                         "name": ".", 
                                         "score": ".", 
                                         "strand": ".",
                                         "seq_id": region_seq_id_list,
                                         "seq_type": region_seq_type_list,
                                         }, 
                                        columns=region_bt_w_anno.column_names,
                                        )
        
        region_bt_w_anno.load_from_dataframe(region_df_w_anno)

        return region_bt_w_anno

    @staticmethod
    def quantify_track_diff(pl_tracks_ref, mn_tracks_ref, 
                            pl_tracks_mut, mn_tracks_mut):
        '''
        Quantify the difference between reference and mutated tracks
        
        Keyword arguments:
        - pl_tracks_ref: plus strand tracks for reference sequences
        - mn_tracks_ref: minus strand tracks for reference sequences
        - pl_tracks_mut: plus strand tracks for mutated sequences
        - mn_tracks_mut: minus strand tracks for mutated sequences

        Returns:
        - track_info_dict: dictionary containing track information
        '''
        track_info_dict = {}

        assert pl_tracks_ref.shape == mn_tracks_ref.shape
        assert pl_tracks_mut.shape == mn_tracks_mut.shape
        assert pl_tracks_ref.shape[0] == 1

        ref_log_total_count = np.log10(pl_tracks_ref.sum() - mn_tracks_ref.sum())
        mut_log_total_count = np.log10(pl_tracks_mut.sum(axis=1) - mn_tracks_mut.sum(axis=1))
        diff_log_total_count = mut_log_total_count - ref_log_total_count

        track_info_dict["ref_log_total_count"] = ref_log_total_count
        track_info_dict["mut_log_total_count"] = mut_log_total_count
        track_info_dict["diff_log_total_count"] = abs(diff_log_total_count)

        return track_info_dict
    
    @staticmethod
    def plot_mutated_vs_ref_predicted_counts(ref_counts, mut_counts, opath=None):
        '''
        Plot the mutated vs. reference predicted counts
        
        Keyword arguments:
        - ref_counts: np.array of reference counts
        - mut_counts: np.array of mutated counts
        - opath: output path for the plot

        Returns:
        - None
        '''
        fig, ax = plt.subplots(1,1, 
                               figsize=(8,8), 
                               )

        x_y_max = max(ref_counts.max(), mut_counts.max())
        ax.plot([0, x_y_max], [0, x_y_max], color="r", linestyle="--")

        ax.scatter(ref_counts, mut_counts, 
                   s=abs(ref_counts - mut_counts)*200,
                   color="k", 
                   )

        ax.set_xlim(0, x_y_max)
        ax.set_ylim(0, x_y_max)
        ax.set_xlabel("Reference Log Total Count")
        ax.set_ylabel("Mutated Log Total Count")
        ax.set_title("Mutated vs. Reference Predicted Counts")

        if opath:
            fig.savefig(opath)

    @staticmethod
    def metaplot_main(args):
        exogeneous_seq = ExogeneousSequences(args.region_file_path, 
                                             args.region_file_type, 
                                             args.fasta,
                                             )
        region_bt = exogeneous_seq.get_region_bed_table()

        region_sizes = region_bt.get_end_locs() - region_bt.get_start_locs()
        elem_size = region_sizes[0]
        num_elem = len(region_bt)

        if not (region_sizes == elem_size).all():
            raise ValueError("Elem sizes should be the same.")

        pl_tracks_arr = np.load(args.signal_tracks_pl)
        mn_tracks_arr = np.load(args.signal_tracks_mn)

        if args.negate_mn_strand:
            mn_tracks_arr = -mn_tracks_arr

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
        exogeneous_seq = ExogeneousSequences(args.region_file_path, 
                                             args.region_file_type, 
                                             args.fasta,
                                             )
        region_bt = exogeneous_seq.get_region_bed_table()

        output_bt = region_bt._clone_empty()
        output_df = pd.DataFrame(columns=output_bt.column_names)

        for region in region_bt.iter_regions():
            seq_name = region["chrom"]

            if re.match(args.regex, seq_name):
                output_df = pd.concat([output_df, 
                                       pd.DataFrame(region.to_dict(), 
                                                    columns=output_df.columns, 
                                                    index=[output_df.shape[0]],
                                                    ), 
                                       ])

        output_bt.write(args.region_out)

        with open(args.fasta, "r") as input_f:
            for record in SeqIO.parse(input_f, "fasta"):
                if re.match(args.regex, record.id):
                    with open(args.fasta_out, "a") as output_f:
                        SeqIO.write(record, output_f, "fasta")

    @staticmethod
    def compare_mutagenesis_main(args):
        exogeneous_seq = ExogeneousSequences(args.region_file_path, 
                                             args.region_file_type, 
                                             args.fasta,
                                             )

        region_bt = exogeneous_seq.get_region_bed_table()
        region_bt_w_anno = ExogeneousTool.match_ref_mut_regex(region_bt,
                                                              args.regex_ref_seqs,
                                                              args.regex_mut_seqs,
                                                              )
        
        pl_tracks = np.load(args.pl_track_npy)
        mn_tracks = np.load(args.mn_track_npy)

        seq_id_list = []
        diff_dict_list = []
        for seq_id in np.unique(region_bt_w_anno.get_region_extra_column("seq_id")):
            ref_logical = (region_bt_w_anno.get_region_extra_column("seq_id") == seq_id) & \
                (region_bt_w_anno.get_region_extra_column("seq_type") == "ref")
            ref_logical = np.array(ref_logical, dtype=bool)
            mut_logical = (region_bt_w_anno.get_region_extra_column("seq_id") == seq_id) & \
                (region_bt_w_anno.get_region_extra_column("seq_type") == "mut")
            mut_logical = np.array(mut_logical, dtype=bool)

            if mut_logical.sum() == 0:
                #TODO: add exception handling here
                continue
            
            track_info_dict = ExogeneousTool.quantify_track_diff(pl_tracks[ref_logical, :], 
                                                                 mn_tracks[ref_logical, :],
                                                                 pl_tracks[mut_logical, :],
                                                                 mn_tracks[mut_logical, :],
                                                                 )
        
            seq_id_list.append(seq_id)
            diff_dict_list.append(track_info_dict)

            with open(os.path.join(args.opath, seq_id + ".json"), "w") as output_f:
                json.dump(track_info_dict, output_f, cls=NumpyEncoder)

        if args.total_count_plot_path:
            ref_counts = np.array([d["ref_log_total_count"] for d in diff_dict_list])

            mut_counts = []
            for d in diff_dict_list:
                largest_diff_ind = np.argmax(d["diff_log_total_count"])
                mut_counts.append(d["mut_log_total_count"][largest_diff_ind])
            
            ExogeneousTool.plot_mutated_vs_ref_predicted_counts(np.array(ref_counts), 
                                                                np.array(mut_counts), 
                                                                opath=args.total_count_plot_path,
                                                                )

    @staticmethod
    def compute_track_correlation_main(args):
        # input check
        if bool(args.mn_track1_npy) ^ bool(args.mn_track2_npy):
            raise ValueError("Either both or none of the minus strand tracks should be provided.")
        mn_track_bool = bool(args.mn_track1_npy)
        
        # load tracks
        track1 = np.load(args.pl_track1_npy)
        track2 = np.load(args.pl_track2_npy)

        if mn_track_bool:
            track1 = np.concatenate([track1, np.load(args.mn_track1_npy)], axis=0)
            track2 = np.concatenate([track2, np.load(args.mn_track2_npy)], axis=0)

        # compute correlation
        if args.jensenshannon:
            js_dist = jensenshannon(track1, track2, axis=1)
            np.save(args.jensenshannon, js_dist)

    @staticmethod
    def main(args):
        if args.subcommand == "metaplot":
            ExogeneousTool.metaplot_main(args)

        elif args.subcommand == "filter":
            ExogeneousTool.filter_main(args)

        elif args.subcommand == "compare_mutagenesis":
            ExogeneousTool.compare_mutagenesis_main(args)

        elif args.subcommand == "compute_track_correlation":
            ExogeneousTool.compute_track_correlation_main(args)

        else:
            raise NotImplementedError


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Exogeneous Tool')
    ExogeneousTool.set_parser(parser)
    args = parser.parse_args()
    sys.exit(ExogeneousTool.main(args))

