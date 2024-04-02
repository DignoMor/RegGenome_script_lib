#!/usr/bin/env python

# overlap TREs and SNPs and return fasta
# and bed files for downstream oligo synthesis

import os
import json
import argparse

import numpy as np
import pandas as pd

from pybedtools import BedTool


def set_parser(parser):
    parser.add_argument("--project_name", 
                        help="Name for the Project.", 
                        dest="project_name", 
                        required=True, 
                        type=str, 
                        )

    parser.add_argument("--tre_inpath", 
                        help="Path to TRE bed file.", 
                        dest="tre_inpath", 
                        required=True, 
                        type=str, 
                        )

    parser.add_argument("--snp_inpath", 
                        help="Path to SNP bed file.", 
                        dest="snp_inpath", 
                        required=True, 
                        type=str, 
                        )

    parser.add_argument("--opath", 
                        help="Output path.", 
                        dest="opath", 
                        required=True, 
                        type=str, 
                        )

    parser.add_argument("--max_length", 
                        help="Maximum Length to be included in final elements.", 
                        dest="max_length", 
                        default=300, 
                        type=int, 
                        )
    
    parser.add_argument("--id_prefix", 
                        help="Prefix for element IDs. Final element ID will be $prefix_$batch_$index.", 
                        dest="id_prefix", 
                        default="TRE", 
                        type=str, 
                        )

    parser.add_argument("--batch_name", 
                        help="Name for element batch. Final element ID will be $prefix_$batch_$index.", 
                        dest="batch_name", 
                        default="E1", 
                        type=str, 
                        )
    
    parser.add_argument("--index_digits", 
                        help="Number of digits for element index.", 
                        dest="index_digits", 
                        default=5, 
                        type=int, 
                        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Overlap TREs and SNPs.")
    set_parser(parser)

    args = parser.parse_args()

    snp_bed = BedTool(args.snp_inpath)
    tre_bed = BedTool(args.tre_inpath)

    tre_intersect_snp_bed = tre_bed.intersect(snp_bed, c=True)
    tre_intersect_snp_df = tre_intersect_snp_bed.to_dataframe()

    # Intersect TREs with risk SNPs
    tre_intersect_snp_df.rename(columns={tre_intersect_snp_df.columns[-1]: "num_snps"}, 
                                inplace=True, 
                                )

    tre_intersect_snp_df['name'] = \
        tre_intersect_snp_df["chrom"] + \
        "_" + tre_intersect_snp_df["start"].agg(str) + \
        "_" + tre_intersect_snp_df["end"].agg(str)
    
    tre_intersect_snp_df["score"] = "."
    tre_intersect_snp_df["strand"] = "."

    tre_w_intersect_df = tre_intersect_snp_df.loc[tre_intersect_snp_df["num_snps"] > 0].copy()

    # Count snp overlap for each candidates elements and summarize
    snp_counts = tre_w_intersect_df["num_snps"].values
    info_dict = dict()

    num_snp_cat, elem_count = np.unique(snp_counts, return_counts=True)
    info_dict["num_snp_histgram"] = {str(n): str(c) for n, c in zip(num_snp_cat, elem_count)} # convert to str for json compatibility
    info_dict["num_snp_elem_list"] = dict()

    for n in num_snp_cat:
        info_dict["num_snp_elem_list"][str(n)] = ",".join(tre_w_intersect_df.loc[tre_w_intersect_df["num_snps"] == n, "name"].values)

    tre_w_intersect_df.drop("num_snps", 
                            axis=1, 
                            inplace=True, 
                            ) # drop the count column to produce standard bed6

    # filter elements by length
    length = tre_w_intersect_df["end"] - tre_w_intersect_df["start"]
    num_tres = length.size
    
    length_filter = length <= args.max_length
    info_dict["length_filter"] = {"passed": str(length_filter.sum()), 
                                  "filtered": str(num_tres - length_filter.sum())}

    filtered_tres_w_intersect_df = tre_w_intersect_df.loc[~length_filter]
    tre_w_intersect_df = tre_w_intersect_df.loc[length_filter]

    # creating index
    tre_id=[str(i).zfill(args.index_digits) for i in tre_w_intersect_df.reset_index(drop=True).index]
    tre_id= [args.id_prefix + "_" + args.batch_name + "_" + s for s in tre_id]

    id2name_table = pd.DataFrame({"id": tre_id, 
                                  "name": tre_w_intersect_df["name"].values, 
                                  }, 
                                 )
    tre_w_intersect_df["name"] = tre_id

    # output
    id2name_table.to_csv(os.path.join(args.opath, args.project_name + ".id_info.tsv",),
                         index=False, 
                         sep="\t", 
                         )
    
    with open(os.path.join(args.opath, args.project_name + ".info.json"), "w") as f:
        f.write(json.dumps(info_dict, indent=4))

    tre_w_intersect_df.to_csv(os.path.join(args.opath, args.project_name + ".tre.bed",),
                              index=False, 
                              header=False, 
                              sep="\t", 
                              )

    filtered_tres_w_intersect_df.to_csv(os.path.join(args.opath, args.project_name + ".tre.filtered.bed",),
                                        index=False, 
                                        header=False, 
                                        sep="\t", 
                                        )
   
