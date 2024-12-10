#!/usr/bin/env python3

import argparse
import sys
import re

from Bio import SeqIO

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
    def metaplot_main(args):
        raise NotImplementedError

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

