#!/usr/bin/env python

import os
import re
import sys
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from umap import UMAP

def set_parser(parser):
    parser.add_argument("--inpath", 
                        "-I", 
                        dest="inpath", 
                        required=True, 
                        help="Path to RSEM result.", 
                        )

    # Output Arguments
    parser.add_argument("--tpm_opath", 
                        "-O", 
                        dest="tpm_opath", 
                        help="Path to TPM table.", 
                        default=None, 
                        )
    
    parser.add_argument("--pca_opath", 
                        dest="pca_opath", 
                        help="Path to PCA plot.", 
                        default=None, 
                        )

    parser.add_argument("--umap_opath", 
                        dest="umap_opath", 
                        help="Path to UMAP plot.", 
                        default=None, 
                        )

    # Optional Arguments
    parser.add_argument("--result_regex", 
                        dest="result_regex", 
                        default="(.*)\.genes\.results", 
                        help="Regular Expression to match for results file, with capture group being the sample name.", 
                        )

    parser.add_argument("--verbose", 
                        dest="verbose", 
                        default=False, 
                        type=lambda s: not s.upper() == "FALSE", 
                        )

def get_file_and_sample_list(inpath, regex_to_match):
    all_file_list = os.listdir(inpath)
    matched_file_list = [os.path.join(inpath, f) for f in all_file_list if re.match(regex_to_match, os.path.basename(f))]
    sample_list = [re.match(regex_to_match, os.path.basename(f))[1] for f in matched_file_list]

    return matched_file_list, sample_list

def get_index_array(file_list, sample_list, verbose=False):
    all_gene_id_list = []
    for file_name in file_list:
        data_df = pd.read_csv(file_name, 
                              sep="\t", 
                              )

        gene_id_list = data_df["gene_id"].values
        all_gene_id_list.append(gene_id_list)
    
    result_list = all_gene_id_list[0]

    for l in all_gene_id_list:
        result_list = np.intersect1d(result_list, l)
    
    if verbose:
        num_gene_id = len(result_list)
        sys.stderr.write("Num Gene IDs: {:d}\n".format(num_gene_id))

        id_num_df = pd.DataFrame({"num_gene_id": [len(l) for l in all_gene_id_list]}, 
                                 index=sample_list, 
                                 )

        sys.stderr.write(id_num_df.to_string())
        sys.stderr.write("\n\n")
    
    return result_list

def make_dim_reduction_plot(data_df, reduction_framework, ax):
    reduction_result = reduction_framework.fit_transform(data_df.values.T)

    for i in range(reduction_result.shape[0]):
        ax.scatter(reduction_result[i, 0],
                   reduction_result[i, 1],  
                   s=10, 
                   color="k", 
                   )
        ax.annotate(data_df.columns[i], 
                    reduction_result[i], 
                    )


def main(args):
    file_list, sample_list = get_file_and_sample_list(args.inpath, args.result_regex)

    index_array = get_index_array(file_list, 
                                  sample_list, 
                                  verbose=args.verbose, 
                                  )

    tpm_df = pd.DataFrame(columns=sample_list, 
                          index=index_array, 
                          )
    
    for file_path, sample_name in zip(file_list, sample_list):
        data_df = pd.read_csv(file_path, 
                              sep="\t", 
                              index_col="gene_id", 
                              )

        tpm_df[sample_name] = data_df.loc[tpm_df.index, "TPM"]

    if args.tpm_opath:
        tpm_df.to_csv(args.tpm_opath)
    
    if args.pca_opath:
        fig, ax = plt.subplots(1,1)
        pca = PCA(n_components=2)
        make_dim_reduction_plot(data_df=tpm_df, 
                                reduction_framework=pca, 
                                ax=ax, 
                                )

        fig.suptitle("PCA")
        fig.tight_layout()
        fig.savefig(args.pca_opath)

    if args.umap_opath:
        fig, ax = plt.subplots(1,1)
        umap = UMAP(n_components=2, 
                    n_neighbors=5, 
                    )
        make_dim_reduction_plot(data_df=tpm_df, 
                                reduction_framework=umap, 
                                ax=ax, 
                                )

        fig.suptitle("UMAP")
        fig.tight_layout()
        fig.savefig(args.umap_opath)
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Summarize RSEM.")

    set_parser(parser=parser)
    args = parser.parse_args()

    main(args)
