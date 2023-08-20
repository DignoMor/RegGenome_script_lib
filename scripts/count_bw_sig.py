#!/usr/bin/env python

import os
import argparse
import numpy as np
import pandas as pd
import pyBigWig

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
                        help="type of regional file. [0]"
                        "(0: TRE-BED file outputed by PINTS)", 
                        default=0, 
                        type=int, 
                        dest="region_file_type", 
                        )

    parser.add_argument("--opath", 
                        help="Output path.", 
                        type=str, 
                        required=True, 
                        dest="opath", 
                        )

def parse_region_input(region_file, file_type):
    '''
    Parse the input file and return a dataframe include region info.
    The returned df could have type-specific additional information that can be 
    used in further analysis, but the followings are always included:
    - "chr"
    - "start_loc"
    - "end_loc"
    
    The start and end loc will following the bed file convention (0-base, half-open).
    The index of the region df will be the unique identifier for the region 
    in the output matrix.
    '''
    if file_type == 0:
        region_df = pd.read_csv(region_file, 
                                sep="\t", 
                                names=["chr", 
                                       "start_loc", 
                                       "end_loc", 
                                       "element_type", 
                                       "fwd_tss",
                                       "rev_tss", 
                                       ], 
                                )

        region_df.index = region_df["chr"] + "_" + region_df["start_loc"].agg(str) + "_" + region_df["end_loc"].agg(str)

    else:
        raise Exception("Unsupported region file type ({}).".format(file_type))
        

    return region_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Extracting Counts from BWs")
    set_parser(parser=parser)
    args = parser.parse_args()

    region_df = parse_region_input(args.region_file_path, 
                                   args.region_file_type, 
                                   )

    count_df = pd.DataFrame(index=region_df.index, 
                            columns=args.sample_names, 
                            )

    

    for sample ,bw_pl_path, bw_mn_path in zip(args.sample_names, args.bw_pls, args.bw_mns):
        bw_pl = pyBigWig.open(bw_pl_path)
        bw_mn = pyBigWig.open(bw_mn_path)

        for region_id, region_info in region_df.iterrows():
            count = np.nan_to_num(bw_pl.values(region_info["chr"], 
                                               region_info["start_loc"], 
                                               region_info["end_loc"], 
                                               ), 
                                  ).sum()
            count -= np.nan_to_num(bw_mn.values(region_info["chr"], 
                                                region_info["start_loc"], 
                                                region_info["end_loc"], 
                                                ), 
                                   ).sum()

            count_df.loc[region_id, sample] = count
        
        bw_pl.close()
        bw_mn.close()

    count_df.to_csv(os.path.join(args.opath, args.job_name + ".count.csv"))
    region_df.to_csv(os.path.join(args.opath, args.job_name + ".region_info.csv"))
