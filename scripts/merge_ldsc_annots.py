#!/usr/bin/env python

import argparse

import pandas as pd

class MergeLDSCAnnots:
    @staticmethod
    def set_parser(parser):
        parser.add_argument('--annot1', 
                            required=True, 
                            help='First annotation file', 
                            )

        parser.add_argument('--annot2',
                            required=True,
                            help='Second annotation file',
                            )
        
        parser.add_argument('--opath',
                            required=True,
                            help='Output file',
                            )

    @staticmethod
    def main(args):
        annot1_df = pd.read_csv(args.annot1, 
                                sep='\t', 
                                compression='gzip',
                                )
        
        annot2_df = pd.read_csv(args.annot2,
                                sep='\t',
                                compression='gzip',
                                )
        
        annot2_df.drop(columns=["CHR", "CM", "BP"],
                       inplace=True,
                       )

        merged_df = pd.merge(annot1_df, 
                             annot2_df, 
                             on='SNP', 
                             how='left', 
                             )
        
        merged_df.to_csv(args.opath,
                         sep='\t',
                         index=False,
                         compression='gzip',
                         )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    MergeLDSCAnnots.set_parser(parser)

    args = parser.parse_args()

    MergeLDSCAnnots.main(args)
