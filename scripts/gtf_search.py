#!/usr/bin/env python

import argparse
import RGTools
import sys

import pandas as pd

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
                        )

    parser.add_argument("--additional_feature_key_value_pair", 
                        help="Specify additional feature to search. "
                             "(eg. gene_type==protein_coding)", 
                        action="append", 
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="Search GTF.")

    set_parser(parser)

    args = parser.parse_args()

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

    start_loc_list = []
    end_loc_list = []
    gene_id_list = []
    chr_name_list = []
    score_list = []
    strand_list = []

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

    bed_df = pd.DataFrame({
        "chr_name": chr_name_list,
        "start": [t - 1 for t in start_loc_list], # bed is 0-based
        "end": [t - 1 for t in end_loc_list],
        "gene_id": gene_id_list,
        "score": score_list,
        "strand": strand_list,
    })

    if args.bed_out == "stdout":
        bed_out = sys.stdout
    else:
        bed_out = args.bed_out

    bed_df.to_csv(bed_out,
                  sep="\t",
                  header=False,
                  index=False,
                  )




