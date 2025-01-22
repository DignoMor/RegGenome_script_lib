#!/usr/bin/env python

import argparse
import sys

import numpy as np

from RGTools.GenomicElements import GenomicElements
from RGTools.BwTrack import BwTrack
from RGTools.utils import str2bool

class GenomicElementTool:
    @staticmethod
    def set_parser(parser):
        subparsers = parser.add_subparsers(dest="subcommand")

        parser_count_bw = subparsers.add_parser("count_bw",
                                                help="Count signal in bigwig files.",
                                                )

        GenomicElementTool.set_parser_count_bw(parser_count_bw)

    @staticmethod
    def set_parser_count_bw(parser):
        GenomicElements.set_parser_genomic_element_region(parser)
        parser.add_argument("--bw_pl",
                            help="Plus strand bigwig file.",
                            required=True,
                            type=str,
                            )

        parser.add_argument("--bw_mn",
                            help="Minus strand bigwig file.",
                            type=str,
                            default=None,
                            )

        parser.add_argument("--single_bw",
                            help="If there is only plus strand bigwig file.",
                            type=str2bool,
                            default=False,
                            )

        parser.add_argument("--ignore_strandness",
                            help="Ignore strandness of the signal.",
                            type=str2bool,
                            default=False,
                            )

        parser.add_argument("--quantification_type",
                            help="Type of quantification.",
                            type=str,
                            default="raw_count",
                            choices=BwTrack.get_supported_quantification_type(),
                            )

        parser.add_argument("--opath",
                            help="Output path for counting.",
                            required=True,
                            type=str,
                            )

    @staticmethod
    def main(args):
        if args.subcommand == "count_bw":
            GenomicElementTool.count_bw_main(args)
        else:
            raise ValueError("Unknown subcommand: {}".format(args.subcommand))

    @staticmethod
    def count_bw_main(args):
        genomic_elements = GenomicElements(region_path=args.region_file_path,
                                           region_file_type=args.region_file_type,
                                           genome_path=args.genome_path,
                                           )
        region_bt = genomic_elements.get_region_bed_table()

        bw_track = BwTrack(bw_pl_path=args.bw_pl,
                           bw_mn_path=args.bw_mn,
                           single_bw=args.single_bw,
                           )

        output_list = []
        for region in region_bt.iter_regions():
            output_list.append(bw_track.count_single_region(region["chrom"],
                                                            region["start"],
                                                            region["end"],
                                                            region["strand"],
                                                            output_type=args.quantification_type,
                                                            ),
                               )
        output_arr = np.array(output_list)

        np.save(args.opath, output_arr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genomic element tool.")
    GenomicElementTool.set_parser(parser)
    args = parser.parse_args()
    GenomicElementTool.main(args)
