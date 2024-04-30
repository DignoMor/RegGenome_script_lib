#!/usr/bin/env python

import argparse
import RGTools
import sys

import pandas as pd

from RGTools.BedTable import BedTable6, BedTable6Plus

def set_parser(parser):
    parser.add_argument("--gtf_path", "-g", 
                        help="path to gtf file.", 
                        )

    parser.add_argument("--bed_out", "-o", 
                        help="Output path.", 
                        )

    parser.add_argument("--general_feature_key_value_pair", 
                        help="Specify general feature to search. "
                             "(eg. chr_name==1; "
                             "Supported General Features: "
                             "chr_name, "
                             "record_source, "
                             "feature_type, "
                             "start_loc, "
                             "end_loc, "
                             "score, "
                             "strand, "
                             "phase"
                             ")", 
                        action="append", 
                        default=[],
                        )

    parser.add_argument("--additional_feature_key_value_pair", 
                        help="Specify additional feature to search. "
                             "(eg. gene_type==protein_coding)", 
                        action="append", 
                        default=[],
                        )
    
    parser.add_argument("--extra_col_general_feature", 
                        help="Extra general feature column to output. "
                             "(eg. feature_type)", 
                        action="append", 
                        default=[],
                        )

    parser.add_argument("--extra_col_additional_feature", 
                        help="Extra additional feature column to output. "
                             "(eg. gene_name)", 
                        action="append", 
                        default=[],
                        )

def key_value_pair_str_parser(key_value_pair_str):
    '''
    Parse key value pair from a string
    
    >>> key_value_pair_str_parser("gene_type==protein_coding")
    ('gene_type', 'protein_coding')
    >>> key_value_pair_str_parser("chr_name==1")
    ('chr_name', '1')
    >>> key_value_pair_str_parser("start_loc==100")
    ('start_loc', 100)
    '''
    key, value = key_value_pair_str.split("==")

    if key in ("start_loc", "end_loc"):
        value = int(value)

    return key, value

def main(args):
    # put searching key-value pairs to dicts
    general_key_value_dict = dict()
    additional_key_value_dict = dict()

    for key_value_str in args.general_feature_key_value_pair:
        search_key, search_value = key_value_pair_str_parser(key_value_str)
        general_key_value_dict[search_key] = search_value

    for key_value_str in args.additional_feature_key_value_pair:
        search_key, search_value = key_value_pair_str_parser(key_value_str)
        additional_key_value_dict[search_key] = search_value

    gtf_file_handle = RGTools.GTFHandle(gtf_path=args.gtf_path)

    # filter the gtf_file
    for search_key, search_value in general_key_value_dict.items():
        gtf_file_handle = gtf_file_handle.filter_by_general_record(search_key, search_value)

    for search_key, search_value in additional_key_value_dict.items():
        gtf_file_handle = gtf_file_handle.filter_by_add_record(search_key, search_value)

    # search the gtf_file
    start_loc_list = []
    end_loc_list = []
    gene_id_list = []
    chr_name_list = []
    score_list = []
    strand_list = []

    extra_col_general_feature_list_dict = {f:[] for f in args.extra_col_general_feature}
    extra_col_add_feature_list_dict = {f:[] for f in args.extra_col_additional_feature}

    for record in gtf_file_handle:
        start_loc = record.search_general_info("start_loc")
        end_loc = record.search_general_info("end_loc")

        gene_id = record.search_add_info("gene_id")
        chr_name = record.search_general_info("chr_name")
        score = record.search_general_info("score")
        strand = record.search_general_info("strand")

        start_loc_list.append(start_loc)
        end_loc_list.append(end_loc)
        gene_id_list.append(gene_id)
        chr_name_list.append(chr_name)
        score_list.append(score)
        strand_list.append(strand)

        for f in extra_col_general_feature_list_dict.keys():
            extra_col_general_feature_list_dict[f].append(record.search_general_info(f))

        for f in extra_col_add_feature_list_dict.keys():
            extra_col_add_feature_list_dict[f].append(record.search_add_info(f))

    bed_df = pd.DataFrame({
        "chr_name": chr_name_list,
        "start": [t - 1 for t in start_loc_list], # GENCODE is 1-based but bed is 0-based
        "end": [t - 1 for t in end_loc_list], # GENCODE is closed but bed is open at the end
        "gene_id": gene_id_list,
        "score": score_list,
        "strand": strand_list,
    })

    for f, v in extra_col_general_feature_list_dict.items():
        bed_df[f] = v
    
    for f, v in extra_col_add_feature_list_dict.items():
        bed_df[f] = v

    if args.bed_out == "stdout":
        bed_out = sys.stdout
    else:
        bed_out = args.bed_out
    
    all_extra_col_names = args.extra_col_general_feature + args.extra_col_additional_feature
    if len(all_extra_col_names) == 0:
        output_bed_table = BedTable6()
    else:
        output_bed_table = BedTable6Plus(extra_column_names=all_extra_col_names, 
                                         extra_column_dtype=["str"] * len(all_extra_col_names),
                                         )
    
    column_map = {"chrom": "chr_name",
                  "name": "gene_id",
                  "score": "score",
                  "strand": "strand",
                  "start": "start",
                  "end": "end",
                  }
    
    for f in args.extra_col_general_feature + args.extra_col_additional_feature:
        column_map[f] = f

    output_bed_table.load_from_dataframe(bed_df, 
                                         column_map=column_map, 
                                         )

    output_bed_table.write(bed_out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="Search GTF.")

    set_parser(parser)

    args = parser.parse_args()

    main(args)
    