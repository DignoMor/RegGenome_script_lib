
import argparse
import sys

sys.path.append("scripts")
sys.path.append(".")
from scripts.GWAS_visualizer import main

if __name__ == "__main__":
    args = argparse.Namespace(gwas_summary_bed_path="sample_data/sample.gwas_summary.bed6",
                              annotation_bed_path=["sample_data/sample.annot.bed3"],
                              annotation_names=["sample_annotation"],
                              plot_chr="chr2",
                              plot_start_loc=127100000,
                              plot_end_loc=127120000,
                              )
    
    main(args)


