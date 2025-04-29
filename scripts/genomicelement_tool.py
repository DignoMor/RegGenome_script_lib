#!/usr/bin/env python

import argparse
import sys

import pandas as pd
import numpy as np

from RGTools.GenomicElements import GenomicElements
from RGTools.exceptions import InvalidBedRegionException
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

        parser_pad_region = subparsers.add_parser("pad_region",
                                                  help="Pad regions. This program differs from "
                                                       "padding_bed.py in that it conserve the " 
                                                       "order of elements in Genomic Elements files.",
                                                  )
        
        GenomicElementTool.set_parser_pad_region(parser_pad_region)

        parser_bed2tssbed = subparsers.add_parser("bed2tssbed",
                                                  help="Convert bed file to TSS bed file.",
                                                  )
        
        GenomicElementTool.set_parser_bed2tssbed(parser_bed2tssbed)

        parser_onehot = subparsers.add_parser("onehot",
                                              help="One-hot encode the sequence. Only support elements of the same size.",
                                              )
        
        GenomicElementTool.set_parser_onehot(parser_onehot)

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

        parser.add_argument("--override_strand",
                            help="overide the strand information in the input file (None if use the input strand info).",
                            type=str,
                            default=None,
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
    def set_parser_pad_region(parser):
        GenomicElements.set_parser_genomic_element_region(parser)
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

        parser.add_argument("--ignore_strand",
                            help="Ignore the strand information in the input file. ",
                            default=False,
                            type=str2bool,
                            )
        
        parser.add_argument("--method_resolving_invalid_region", 
                            help="Method to resolve invalid region after padding.", 
                            type=str,
                            default="fallback",
                            choices=["raise", "fallback", "drop"],
                            )

    @staticmethod
    def set_parser_bed2tssbed(parser):
        GenomicElements.set_parser_genomic_element_region(parser)
        parser.add_argument("--opath",
                            help="Output path for the TSS BED file",
                            type=str,
                            required=True,
                            )

        parser.add_argument("--output_site", 
                            help="The site for output. [TSS] ({})".format(
                                ", ".join(GenomicElementTool.get_bed2tssbed_output_site_types()), 
                            ), 
                            default="TSS",
                            type=str, 
                            )

    def set_parser_onehot(parser):
        GenomicElements.set_parser_genome(parser)
        GenomicElements.set_parser_genomic_element_region(parser)

        parser.add_argument("--opath",
                            help="Output path for the one-hot encoded sequence.",
                            type=str,
                            required=True,
                            )


    @staticmethod
    def get_bed2tssbed_output_site_types():
        return ["TSS", "center"]

    @staticmethod
    def get_bed2tssbed_site_coord(region, output_site_type):
        if output_site_type == "TSS":
            if region["strand"] == "+":
                site = region["start"]
            elif region["strand"] == "-":
                site = region["end"] - 1
        elif output_site_type == "center":
            site = int((region["start"] + region["end"]) // 2)
        
        return site

    @staticmethod
    def StrandInputType(string):
        '''
        Input type for strand information.
        Converts "None" to None.

        Parameters:
        - string: str input for strand info.
        '''
        if string == "None":
            return None
        else:
            return string

    @staticmethod
    def main(args):
        if args.subcommand == "count_bw":
            GenomicElementTool.count_bw_main(args)
        elif args.subcommand == "pad_region":
            GenomicElementTool.pad_region_main(args)
        elif args.subcommand == "bed2tssbed":
            GenomicElementTool.bed2tssbed_main(args)
        elif args.subcommand == "onehot":
            GenomicElementTool.onehot_main(args)
        else:
            raise ValueError("Unknown subcommand: {}".format(args.subcommand))

    @staticmethod
    def pad_region_main(args):
        genomic_elements = GenomicElements(region_path=args.region_file_path,
                                           region_file_type=args.region_file_type,
                                           fasta_path=None, 
                                           )
        
        region_bt = genomic_elements.get_region_bed_table()

        output_dict_list = []

        for i, region in enumerate(region_bt.iter_regions()):
            try:
                new_region = region.pad_region(upstream_padding=args.upstream_pad,
                                               downstream_padding=args.downstream_pad,
                                               ignore_strand=args.ignore_strand,
                                               )
            except InvalidBedRegionException as e:
                if args.method_resolving_invalid_region == "raise":
                    raise e
                elif args.method_resolving_invalid_region == "fallback":
                    new_region = region
                elif args.method_resolving_invalid_region == "drop":
                    continue
                else:
                    raise ValueError(f"Unknown method to resolve invalid region: {args.method_resolving_invalid_region}")

            output_dict_list.append(new_region.to_dict())
        
        output_region_bt = region_bt._clone_empty()
        output_region_bt.load_from_dataframe(pd.DataFrame(output_dict_list,
                                                          columns=region_bt.column_names,
                                                          ))
            
        output_region_bt.write(args.opath)

    @staticmethod
    def count_bw_main(args):
        genomic_elements = GenomicElements(region_path=args.region_file_path,
                                           region_file_type=args.region_file_type,
                                           fasta_path=None, 
                                           )
        region_bt = genomic_elements.get_region_bed_table()

        bw_track = BwTrack(bw_pl_path=args.bw_pl,
                           bw_mn_path=args.bw_mn,
                           single_bw=args.single_bw,
                           )

        output_list = []
        for region in region_bt.iter_regions():
            if args.override_strand:
                strand = args.override_strand
            else:
                strand = "." if args.region_file_type == "bed3" else region["strand"]

            output_list.append(bw_track.count_single_region(region["chrom"],
                                                            region["start"],
                                                            region["end"],
                                                            strand, 
                                                            output_type=args.quantification_type,
                                                            ),
                               )

        output_arr = np.array(output_list)

        np.save(args.opath, output_arr)
    
    @staticmethod
    def bed2tssbed_main(args):
        genomic_elements = GenomicElements(region_path=args.region_file_path,
                                           region_file_type=args.region_file_type,
                                           fasta_path=None, 
                                           )
        region_bt = genomic_elements.get_region_bed_table()

        output_region_list = []
        for region in region_bt.iter_regions():
            output_coord = GenomicElementTool.get_bed2tssbed_site_coord(region, args.output_site)
            
            output_region_dict = region.to_dict()
            output_region_dict["start"] = output_coord
            output_region_dict["end"] = output_coord + 1
            output_region_list.append(output_region_dict)

        output_df = pd.DataFrame(output_region_list, 
                                 columns=region_bt.column_names, 
                                 )
        
        output_bt = region_bt._clone_empty()
        output_bt.load_from_dataframe(output_df)
        output_bt.write(args.opath)

    @staticmethod
    def onehot_main(args):
        genomic_elements = GenomicElements(region_path=args.region_file_path,
                                           region_file_type=args.region_file_type,
                                           fasta_path=args.fasta_path, 
                                           )

        output_arr = genomic_elements.get_all_region_one_hot()
        np.save(args.opath, output_arr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="WARNING: DEPRECATED! " 
                                                 "Please use https://github.com/DignoMor/GenomicElementTool "
                                                 "for the same functionality.", 
                                    )
    GenomicElementTool.set_parser(parser)
    args = parser.parse_args()
    GenomicElementTool.main(args)
