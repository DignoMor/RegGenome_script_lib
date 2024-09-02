#!/usr/bin/env python
import argparse
import pyBigWig

import pandas as pd
import numpy as np

from RGTools.BedTable import BedTable6

class IdentifyPausingRegion:
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--tss_bed", 
                            help="Input path for TSS reference file.", 
                            type=str, 
                            required=True, 
                            )
        
        parser.add_argument("--l_pad", 
                            help="Left padding for pausing region. "
                                 "Padded region will be the search space for pausing region.", 
                            type=int, 
                            default=500,
                            )

        parser.add_argument("--r_pad", 
                            help="Right padding for pausing region. "
                                 "Padded region will be the search space for pausing region.", 
                            type=int, 
                            default=499,
                            )
        
        parser.add_argument("--target_size", 
                            help="Target size for pausing region. "
                                 "Pausing region will be the target size region with the highest signal density.", 
                            type=int, 
                            default=250,
                            )
        
        parser.add_argument("--search_step_size",
                            help="Step size for searching pausing region.",
                            type=int,
                            default=10,
                            )

        parser.add_argument("--bw_pl", 
                            help="Input path for plus strand bigwig file for PROseq assay.", 
                            type=str,
                            required=True,
                            )

        parser.add_argument("--bw_mn", 
                            help="Input path for minus strand bigwig file for PROseq assay.", 
                            type=str,
                            required=True,
                            )
        
        parser.add_argument("--opath",
                            help="Output path for pausing region bed file.",
                            type=str,
                            required=True,
                            )

    @staticmethod
    def max_subset_sig(sig: np.ndarray, target_size: int, 
                       step_size: int=10):
        '''
        Find the subset of the signal with the highest signal density.
        euqals are broken by using the leftmost subset.

        Return:
        - start position (in array index) of the subset.
        - end position (in array index) of the subset.
        - read per kilobase (rpk) of the subset.

        Keyword arguments:
        - sig: signal array.
        - target_size: target size of the subset.
        '''
        max_rpk = 0
        max_start = 0
        max_end = max_start + target_size

        for i in range(0, len(sig) - target_size, step_size):
            sub_sig = sig[i:i+target_size]
            rpk = np.sum(sub_sig) / target_size * 1e3
            if rpk > max_rpk:
                max_rpk = rpk
                max_start = i
                max_end = i + target_size

        return (max_start, max_end, max_rpk)

    @staticmethod
    def search_region(bw_pl: "pyBigWig.bigWigFile",
                      bw_mn: "pyBigWig.bigWigFile",
                      chrom: str, 
                      start: int, 
                      end: int, 
                      strand: str, 
                      target_size: int,
                      search_step_size=10, 
                      ):
        '''
        Search pausing region in the given region.
        The region definition follows the bed convention.

        Return: 
        - start position of the pausing region.
        - end position of the pausing region.
        - read per kilobase (rpk) of the pausing region.

        Keyword arguments:
        - bw_pl: pyBigWig.bigWigFile object for plus strand signal.
        - bw_mn: pyBigWig.bigWigFile object for minus strand signal.
        - chrom: chromosome name.
        - start: start position of the region, inclusive.
        - end: end position of the region, non-inclusive.
        - strand: strand of the gene so that search focus on the particular strand.
        - target_size: target size of the pausing region.
        - search_step_size: step size for searching pausing region.
        '''
        if strand == "+":
            bw = bw_pl
            sig = bw.values(chrom, start, end)
            sig = np.nan_to_num(sig)
        elif strand == "-":
            bw = bw_mn
            sig = bw.values(chrom, start, end)
            sig = -np.nan_to_num(sig)
            sig = np.flip(sig)
        else:
            raise ValueError("Invalid strand")

        start2upstream, end2upstream, rpk = IdentifyPausingRegion.max_subset_sig(sig,
                                                                          target_size=target_size,
                                                                          step_size=search_step_size,
                                                                          )
        
        if strand == "+":
            start = start + start2upstream
            end = start + target_size
        elif strand == "-":
            end = end - start2upstream
            start = end - target_size
        else:
            raise ValueError("Invalid strand")

        return (start, end, rpk)

    @staticmethod
    def main(args):
        tss_bt = BedTable6()
        tss_bt.load_from_file(args.tss_bed)

        bw_pl = pyBigWig.open(args.bw_pl)
        bw_mn = pyBigWig.open(args.bw_mn)

        output_df = pd.DataFrame(columns=["chrom", "start", "end", "name", "score", "strand"])

        for region in tss_bt.iter_regions():
            search_start = region["start"] - args.l_pad
            search_end = region["end"] + args.r_pad + 1
            result_start, result_end, rpk = IdentifyPausingRegion.search_region(bw_pl, 
                                                                                bw_mn, 
                                                                                chrom=region["chrom"], 
                                                                                start=search_start, 
                                                                                end=search_end, 
                                                                                strand=region["strand"],
                                                                                target_size=args.target_size,
                                                                                search_step_size=args.search_step_size,
                                                                                )
            
            output_df.loc[len(output_df)] = {
                "chrom": region["chrom"],
                "start": result_start,
                "end": result_end,
                "name": region["name"],
                "score": rpk,
                "strand": region["strand"],
            }

        output_bt = BedTable6()
        output_bt.load_from_dataframe(output_df)
        output_bt.write(args.opath)

        bw_pl.close()
        bw_mn.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify pausing region")

    IdentifyPausingRegion.set_parser(parser)

    args = parser.parse_args()

    IdentifyPausingRegion.main(args)
