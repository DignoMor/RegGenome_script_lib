#!/usr/bin/env python

import os
import argparse
import numpy as np
import pandas as pd
import pyBigWig

from RGTools.utils import str2bool
from RGTools.BedTable import BedTable3, BedTable6, BedTable6Plus

#TODO: add BedTableTRE support
REGION_FILE_SUFFIX2CLASS_DICT = {
    "bed3": BedTable3,
    "bed6": BedTable6,
    "bed6plus": BedTable6Plus, 
}

def set_parser(parser):
    parser.add_argument("--job_name", 
                        help="Name of the job (also output file prefix).", 
                        required=True,
                        type=str, 
                        dest="job_name", 
                        )

    parser.add_argument("--sample", 
                        help="Sample name.", 
                        action="append", 
                        required=True, 
                        type=str, 
                        dest="sample_names", 
                        )

    parser.add_argument("--bw_pl", 
                        help="bigwig file for plus strand (must match sample order).", 
                        action="append", 
                        required=True, 
                        type=str, 
                        dest="bw_pls", 
                        )

    parser.add_argument("--bw_mn", 
                        help="bigwig file for minus strand (must match sample order).", 
                        action="append", 
                        required=True, 
                        type=str, 
                        dest="bw_mns", 
                        )
    
    parser.add_argument("--region_file_path", 
                        help="Path to the region file.", 
                        type=str, 
                        required=True, 
                        dest="region_file_path", 
                        )
    
    parser.add_argument("--region_file_type.", 
                        help="type of regional file. [bed3]"
                        "(options: {})".format(", ".join(REGION_FILE_SUFFIX2CLASS_DICT.keys())), 
                        default=0, 
                        type=str, 
                        dest="region_file_type", 
                        )

    parser.add_argument("--opath", 
                        help="Output path.", 
                        type=str, 
                        required=True, 
                        dest="opath", 
                        )
    
    parser.add_argument("--ignore_strandness", 
                        help="If to ignore strandness.", 
                        type=str2bool, 
                        default=False, 
                        dest="ignore_strandness", 
                        )

def parse_region_input(region_file, file_type):
    '''
    Parse the input file and return a dataframe include region info.
    The returned df could have type-specific additional information that can be 
    used in further analysis, but the followings are always included:
    - "chrom"
    - "start"
    - "end"
    - "name"
    - "score"
    - "strand"
    
    The start and end loc will following the bed file convention (0-base, half-open).
    The index of the region df will be the unique identifier for the region 
    in the output matrix.
    '''
    if not file_type in REGION_FILE_SUFFIX2CLASS_DICT.keys():
        raise Exception("Unsupported region file type ({}).".format(file_type))
    
    region_bed_table = REGION_FILE_SUFFIX2CLASS_DICT[file_type]()
    region_bed_table.load_from_file(region_file)
    
    if file_type == "bed3":
        new_region_bed_table = BedTable6()
        new_region_bed_table.load_from_BedTable3(region_bed_table)
        region_bed_table = new_region_bed_table

    #TODO: Switch from df to BedTables
    region_df = region_bed_table.to_dataframe()
    region_df.fillna(".", inplace=True)
    region_df.index = region_df["chrom"] + "_" + region_df["start"].transform(str) + "_" + region_df["end"].transform(str)

    return region_df

def main(args):
    region_df = parse_region_input(args.region_file_path, 
                                   args.region_file_type, 
                                   )

    if args.ignore_strandness:
        region_df["strand"] = "."

    count_df = pd.DataFrame(index=region_df.index, 
                            columns=args.sample_names, 
                            dtype=int, 
                            )


    for sample ,bw_pl_path, bw_mn_path in zip(args.sample_names, args.bw_pls, args.bw_mns):
        bw_pl = pyBigWig.open(bw_pl_path)
        bw_mn = pyBigWig.open(bw_mn_path)

        for region_id, region_info in region_df.iterrows():
            pl_count = np.nan_to_num(bw_pl.values(region_info["chrom"], 
                                                  region_info["start"], 
                                                  region_info["end"], 
                                                  ), 
                                     ).sum()
            mn_count = -np.nan_to_num(bw_mn.values(region_info["chrom"], 
                                                   region_info["start"], 
                                                   region_info["end"], 
                                                   ), 
                                      ).sum()

            # Strandness should be encoded in the input 
            # Otherwise both strands are counted and summed
            # TODO: add to documentation
            if region_info["strand"] == "+":
                count = pl_count
            elif region_info["strand"] == "-":
                count = mn_count
            elif region_info["strand"] == ".":
                count = pl_count + mn_count
            else:
                raise Exception("Unexpected value for strand field ({})".format(region_info["strand"]))

            count_df.loc[region_id, sample] = count
        
        bw_pl.close()
        bw_mn.close()

    count_df.to_csv(os.path.join(args.opath, args.job_name + ".count.csv"))
    region_df.to_csv(os.path.join(args.opath, args.job_name + ".region_info.csv"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Extracting Counts from BWs")
    set_parser(parser=parser)
    args = parser.parse_args()

    main(args)
