
import argparse

import pandas as pd

class AnnotateFragments:
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--fragments", "-F",
                            help="Fragments file.",
                            required=True,
                            )
        
        parser.add_argument("--annotation", 
                            help="Annotation file.",
                            required=True,
                            )
        
        parser.add_argument("--oheader", "-O",
                            help="Output header.",
                            required=True,
                            )
        
    @staticmethod
    def main(args):
        fragment_df = pd.read_csv(args.fragments, 
                                  sep="\t", 
                                  header=None, 
                                  names=["chrom", "start", "end", "cell_id", "count"],
                                  index_col="cell_id",
                                  )
        anno_df = pd.read_csv(args.annotation, 
                              sep="\t", 
                              index_col=0, 
                              header=None, 
                              names=["cell_id", "MajorCellType"],
                              )

        uniq_annos = anno_df["MajorCellType"].unique()
        cell_id_set_by_anno = [set(anno_df[anno_df["MajorCellType"] == anno].index) for anno in uniq_annos]

        for anno, cell_id_set in zip(uniq_annos, cell_id_set_by_anno):
            fragment_df_sub = fragment_df.loc[fragment_df.index.isin(cell_id_set)]
            fragment_df_sub.reset_index(inplace=True)
            fragment_df_sub[["chrom", "start", "end", "cell_id", "count"]].to_csv(args.oheader + "." + anno + ".tsv", 
                                                                                  sep="\t", 
                                                                                  index=False, 
                                                                                  header=False, 
                                                                                  )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    AnnotateFragments.set_parser(parser)
    args = parser.parse_args()
    AnnotateFragments.main(args)
        
