
import argparse
import shutil
import sys
import os

sys.path.append("scripts")
sys.path.append(".")
from scripts.gtf_search import main as gtf_search_main
from scripts.GWAS_visualizer import main

if __name__ == "__main__":
    temp_dir = "temp"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gtf_search_args = argparse.Namespace(gtf_path="sample_data/sample.GENCODE.gtf",
                                         bed_out=os.path.join(temp_dir, "sample.transcript_w_gene_name.bed"),
                                         general_feature_key_value_pair=["feature_type==transcript"], 
                                         additional_feature_key_value_pair=["gene_name==BIN1", "gene_type==protein_coding"], 
                                         extra_col_additional_feature=["gene_name"],
                                         extra_col_general_feature=[],
                                         )
    gtf_search_main(gtf_search_args)

    args = argparse.Namespace(gwas_summary_bed_path="sample_data/sample.gwas_summary.bed6",
                              annotation_bed_path=["sample_data/sample1.annot.bed3", 
                                                   os.path.join(temp_dir, "sample.transcript_w_gene_name.bed"), 
                                                   ],
                              annotation_names=["sample_annotation", 
                                                "GENCODE_transcripts", 
                                                ],
                              annotation_file_type=["bed3",
                                                    "bed6Plus_with_gene_name",
                                                    ],
                              plot_chr="chr2",
                              plot_start_loc=127000000,
                              plot_end_loc=127130000,
                              )
    
    main(args)

    shutil.rmtree(temp_dir)
