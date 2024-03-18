#!/usr/bin/env python

import argparse
import re
import pandas as pd

import RGTools

def set_parser(parser):
    parser.add_argument("--gtf_path", "-g",
                        help="Path to gtf file.",
                        )

    parser.add_argument("--bed_out", '-o',
                        help="output path for bed results.",
                        )

    parser.add_argument("--window_size", "-w",
                        help="The window size [250-250]. The default is to return 501 bp window with tss in middle.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Gen TSS bed from GTF.")
    set_parser(parser=parser)
    args = parser.parse_args()

    input_gtf = RGTools.GTFHandle(args.gtf_path)
    protein_coding_gtf = input_gtf.filter_by_add_record("gene_type", "protein_coding")
    protein_coding_transcripts_gtf = protein_coding_gtf.filter_by_general_record("feature_type", "transcript")

    tss_list = []
    tss_plus_one_list = []
    gene_id_list = []
    chr_name_list = []
    score_list = []
    strand_list = []

    for record in protein_coding_transcripts_gtf:
        if record.search_general_info("strand") == "+":
            tss = record.search_general_info("start_loc")
        else:
            tss = record.search_general_info("end_loc") - 1 # GTF is closed at start and open at end

        gene_id = record.search_add_info("gene_id")
        chr_name = record.search_general_info("chr_name")
        score = record.search_general_info("score")
        strand = record.search_general_info("strand")

        tss_list.append(tss)
        gene_id_list.append(gene_id)
        chr_name_list.append(chr_name)
        score_list.append(score)
        strand_list.append(strand)

    match = re.search("([0-9]*)-([0-9]*)", args.window_size)
    if not match:
        raise Exception("Unrecognized window size: {}".format(args.window_size))

    l_pad = int(match.group(1))
    r_pad = int(match.group(2))

    bed_df = pd.DataFrame({
        "chr_name": chr_name_list,
        "start": [t - l_pad - 1 for t in tss_list], # bed is 0-based
        "end": [t + r_pad for t in tss_list],
        "gene_id": gene_id_list,
        "score": score_list,
        "strand": strand_list,
    })

    bed_df.to_csv(args.bed_out,
                  sep="\t",
                  header=False,
                  index=False,
                  )

