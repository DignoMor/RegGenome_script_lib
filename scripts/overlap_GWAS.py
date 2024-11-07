#!/usr/bin/env python

# Overlap GWAS with genomic annotations and return overlaps 
# of different tiers

import argparse
import os

import numpy as np

from RGTools.BedTable import BedTable3, BedTable6

class OverlapGWAS:
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--job_name",
                            help="Job name.",
                            dest="job_name",
                            required=True,
                            type=str,
                            )

        parser.add_argument("--region_path", 
                            help="Path to region bed file.", 
                            dest="region_path", 
                            required=True, 
                            type=str, 
                            )
        
        parser.add_argument("--gwas_path",
                            help="Path to GWAS bed file.",
                            dest="gwas_path",
                            required=True,
                            type=str,
                            )
        
        parser.add_argument("--opath",
                            help="Output path.",
                            dest="opath",
                            required=True,
                            type=str,
                            )
        
        parser.add_argument("--tier_score_cutoff",
                            help="Tier score cutoff.",
                            dest="tier_score_cutoff",
                            type=float,
                            action="append", 
                            )
        
        parser.add_argument("--region_file_type", 
                            help="File type of region file [bed6]. (Opations: {})".format(OverlapGWAS.get_region_file_types()),
                            dest="region_file_type",
                            default="bed6", 
                            choices=OverlapGWAS.get_region_file_types(),
                            type=str, 
                            )
        
        

    @staticmethod
    def get_region_file_types():
        return ["bed6", "bed3"]

    @staticmethod
    def load_region_file(region_path, region_file_type):
        '''
        Load region file from path
        
        Keyword arguments:
        - region_path: path to region file
        - region_file_type: file type of region file

        Returns:
        - BedTable object
        '''
        if region_file_type == "bed6":
            bt = BedTable6()
        elif region_file_type == "bed3":
            bt = BedTable3()

        bt.load_from_file(region_path)

        return bt

    @staticmethod
    def main(args):

        region_bt = OverlapGWAS.load_region_file(args.region_path, args.region_file_type)
        gwas_bt = OverlapGWAS.load_region_file(args.gwas_path, "bed6")

        region_hit_index_by_tier = []

        for tier, score_cutoff in enumerate(args.tier_score_cutoff):
            tier_gwas_bt = gwas_bt.apply_logical_filter(gwas_bt.get_region_scores() < score_cutoff)
            hit_region_ind = set()

            for region in tier_gwas_bt.iter_regions():
                hit_region_ind.update(region_bt.search_region(region["chrom"], 
                                                              region["start"], 
                                                              region["end"], 
                                                              overlapping_base=1, 
                                                              ))
            
            if not len(region_hit_index_by_tier) == 0:
                for ind in np.concatenate(region_hit_index_by_tier):
                    hit_region_ind.discard(ind)
            
            region_hit_index_by_tier.append(list(hit_region_ind))

        for tier, hit_index in enumerate(region_hit_index_by_tier):
            tier_hit_bt = region_bt.apply_logical_filter(hit_index)
            tier_hit_bt.write(os.path.join(args.opath, "{}.tier{:d}.bed".format(args.job_name, tier+1)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Overlap GWAS with genomic annotations and return overlaps of different tiers")
    OverlapGWAS.set_parser(parser)
    args = parser.parse_args()
    OverlapGWAS.main(args)
