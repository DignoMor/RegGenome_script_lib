#!/usr/bin/env python

# convert RSEM output to count table
#TODO: incorporate this script to count table tool and remove this script

import os
import argparse

import pandas as pd
import numpy as np

class RSEM2CountTable:
    def __init__(self):
        pass

    @staticmethod
    def set_parser(parser):
        parser.add_argument("--sample",
                            help="Name of the sample.", 
                            dest="sample", 
                            action="append", 
                            required=True, 
                            )

        parser.add_argument("--rsem_output", 
                            help="Path to the RSEM output.", 
                            dest="rsem_output", 
                            action="append", 
                            required=True, 
                            )
        
        parser.add_argument("--gene_id_list",
                            help="A id list file containing gene ids.", 
                            dest="gene_id_list", 
                            required=True, 
                            )
        
        parser.add_argument("--opath", "-O", 
                            help="Path to the output count table.", 
                            dest="opath", 
                            required=True, 
                            )

        parser.add_argument("--count_type",
                            help="Type of the count, either 'expected_count' or 'TPM'.", 
                            dest="count_type", 
                            default="expected_count", 
                            )
    
    @staticmethod
    def read_rsem_output(rsem_output):
        rsem_table = pd.read_csv(rsem_output,
                                 sep="\t",
                                 index_col=0,
                                 )
        return rsem_table

    @staticmethod
    def read_list(gene_id_list):
        with open(gene_id_list, "r") as f:
            gene_id_list = f.read().splitlines()
        return gene_id_list
    
    @staticmethod
    def get_count_by_id(rsem_table, gene_id_list, count_type):
        counts = []
        for gene_id in gene_id_list:
            if gene_id in rsem_table.index:
                count = rsem_table.loc[gene_id, count_type]
            else:
                count = np.nan
            counts.append(count)
        return np.array(counts)

    @staticmethod
    def save_count_table(count_table, opath):
        count_table.to_csv(opath)
        return None

    @staticmethod
    def main(args):
        rsem_tables = [RSEM2CountTable.read_rsem_output(rsem_output) for rsem_output in args.rsem_output]
        gene_ids = RSEM2CountTable.read_list(args.gene_id_list)

        count_dict = {s: RSEM2CountTable.get_count_by_id(rsem_table, gene_ids, args.count_type) for s, rsem_table in zip(args.sample, rsem_tables)}
        count_table = pd.DataFrame(count_dict, index=gene_ids)

        RSEM2CountTable.save_count_table(count_table, args.opath)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Convert RSEM output to count table.")
    RSEM2CountTable.set_parser(parser)
    args = parser.parse_args()
    RSEM2CountTable.main(args)
