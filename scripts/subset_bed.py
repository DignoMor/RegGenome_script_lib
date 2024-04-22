#!/usr/bin/env python

import argparse

from RGTools.BedTable import BedTable3, BedTable6

def set_parser(parser:argparse.ArgumentParser):
    parser.add_argument("--input_bed_path",
                        help="path to input bed file.",
                        )
    
    parser.add_argument("--output_bed_path",
                        help="path to output bed file.",
                        )
    
    parser.add_argument("--chr_name",
                        help="Chromosome name to subset.",
                        )

    parser.add_argument("--start_loc",
                        help="Start location to subset.",
                        type=int,
                        )
    
    parser.add_argument("--end_loc",
                        help="End location to subset.",
                        type=int,
                        )

    parser.add_argument("--file_type",
                        help="File type of input bed file.",
                        choices=["bed3", "bed6"],
                        )
    
def main(args):
    if args.file_type == "bed3":
        bed_table = BedTable3()
    elif args.file_type == "bed6":
        bed_table = BedTable6()
    else:
        raise ValueError("Invalid file type.")
    
    bed_table.load_from_file(args.input_bed_path)
    bed_table.region_subset(args.chr_name, args.start_loc, args.end_loc).write(args.output_bed_path)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog="Subset bed.")

    set_parser(parser)

    args = parser.parse_args()

    main(args)
