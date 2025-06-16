#! /usr/bin/env python

import argparse

import pandas as pd
import numpy as np

class SplitFragments:
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--fragments", 
                            type=str, 
                            required=True, 
                            help="Fragments file", 
                            )

        parser.add_argument("--label_file", 
                            type=str, 
                            required=True, 
                            help="Label file.", 
                            )

        parser.add_argument("--oheader", 
                            type=str, 
                            required=True, 
                            help="Output file header.", 
                            )

    @staticmethod
    def main(args):
        if args.fragments.endswith(".gz"):
            fragment_df = pd.read_csv(args.fragments, 
                                      sep="\t", 
                                      names=["chrom", "start", "end", "barcode", "support"], 
                                      compression="gzip", 
                                      )
        else:
            fragment_df = pd.read_csv(args.fragments, 
                                    sep="\t", 
                                    names=["chrom", "start", "end", "barcode", "support"], 
                                    )

        label_df = pd.read_csv(args.label_file, 
                               sep="\t", 
                               names=["cell_type"], 
                               index_col=0, 
                               )
        
        label_dict = {k:v for k, v in zip(label_df.index, label_df["cell_type"])}

        for cell_type in label_df.loc[:, "cell_type"].unique():
            cell_type_logical_list = []
            for barcode in fragment_df["barcode"]:
                if barcode in label_dict.keys():
                    cell_type_logical_list.append(label_dict[barcode] == cell_type)
                else:
                    cell_type_logical_list.append(False)

            cell_type_df = fragment_df.loc[np.array(cell_type_logical_list)]
            cell_type_df.to_csv(args.oheader + "." + cell_type + ".fragments.tsv", 
                                sep="\t", 
                                index=False, 
                                header=False,
                                )
                


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    SplitFragments.set_parser(parser)
    args = parser.parse_args()
    SplitFragments.main(args)
