#!/usr/bin/env python

# Lift over SNP bed from one genome version to another
# Uses Ensembl REST API to get the rsid of each SNP, 
# and then map rsid to target genome version


import argparse

import pandas as pd

from RGTools.SNP_utils import EnsemblRestSearch
from RGTools.BedTable import BedTable6Plus
from RGTools.logging import Logger


class LiftOverSNPBed:
    @staticmethod
    def SNPBedTable():
        bt = BedTable6Plus(extra_column_names=["bases"], 
                           extra_column_dtype=[str], 
                           )
        return bt
    
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--bed_in", "-I",
                            help="Path to input SNP bed file.",
                            required=True, 
                            )
        
        parser.add_argument("--bed_out", '-O',
                            help="output path for bed results.",
                            required=True, 
                            )
        
        parser.add_argument("--from_genome",
                            help="The genome version of the input bed file. "
                                 "[hg19] ({})".format(", ".join(LiftOverSNPBed.get_supported_genomes())),
                            default="hg19", 
                            )
        
        parser.add_argument("--to_genome",
                            help="The genome version to lift over to. "
                                 "[hg38] ({})".format(", ".join(LiftOverSNPBed.get_supported_genomes())),
                            required="hg38", 
                            )
        
        parser.add_argument("--log_file", "-L",
                            help="Path to log file.",
                            default="stderr", 
                            )

    @staticmethod
    def get_supported_genomes():
        sample_rest_search = EnsemblRestSearch()
        return list(sample_rest_search.genome_version2url_dict.keys())
    
    @staticmethod
    def main(args):
        logger = Logger(args.log_file)

        input_bt = LiftOverSNPBed.SNPBedTable()
        input_bt.load_from_file(args.bed_in)

        ensembl_rest_source_assembly = EnsemblRestSearch(genome_version=args.from_genome)
        ensembl_rest_target_assembly = EnsemblRestSearch(genome_version=args.to_genome)


        output_chrom_list = []
        output_start_list = []
        output_end_list = []
        output_name_list = []
        output_score_list = []
        output_strand_list = []
        output_bases_list = []

        for region in input_bt.iter_regions():
            try:
                rsids = ensembl_rest_source_assembly.get_rsid_from_location(region["chrom"], 
                                                                            region["start"], 
                                                                            ) 
                rsid, target_region_simple_info = ensembl_rest_target_assembly.prioritize_rsids(rsids)
            except:
                logger.take_log("Failed to get rsid for region: {}\t{:d}\t{:d}.".format(region["chrom"],
                                                                                        region["start"],
                                                                                        region["end"],
                                                                                        ))
                continue

            if rsid is None:
                logger.take_log("Failed to get rsid for region: {}\t{:d}\t{:d}.".format(region["chrom"],
                                                                                        region["start"],
                                                                                        region["end"],
                                                                                        ))
                continue

            output_chrom_list.append(target_region_simple_info["chrom"])
            output_start_list.append(target_region_simple_info["start"])
            output_end_list.append(target_region_simple_info["end"])
            output_name_list.append(rsid)
            output_score_list.append(region["score"])
            output_strand_list.append(region["strand"])
            output_bases_list.append(region["bases"])

        output_bt = LiftOverSNPBed.SNPBedTable()
        output_df = pd.DataFrame({"chrom": output_chrom_list,
                                  "start": output_start_list,
                                  "end": output_end_list,
                                  "name": output_name_list,
                                  "score": output_score_list,
                                  "strand": output_strand_list,
                                  "bases": output_bases_list,
                                  })
        output_bt.load_from_dataframe(output_df)

        output_bt.write(args.bed_out)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    LiftOverSNPBed.set_parser(parser)
    args = parser.parse_args()
    LiftOverSNPBed.main(args)
