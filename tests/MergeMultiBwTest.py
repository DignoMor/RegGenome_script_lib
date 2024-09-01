
import os
import sys
import shutil
import pyBigWig
import argparse
import unittest

import pandas as pd

sys.path.append("scripts")
from scripts.merge_multi_bw import MergeMultiBw
from scripts.merge_bw import merge_bw_files
#TODO: implement tests independent of the legacy script

class MergeMultiBwTest(unittest.TestCase):
    def setUp(self) -> None:
        self.__test_dir = "merge_multi_bw_test"
        self.__chrom_size_path = os.path.join(self.__test_dir, "chrom.size")

        self.chrom_size_df = pd.DataFrame({"chrom": ["chr3", "chr4", "chrFake"],
                                           "size": [198295559, 190214555, 1000000],
                                           },
                                          columns=["chrom", "size"],
                                          )

        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        self.chrom_size_df.to_csv(self.__chrom_size_path,
                                  sep="\t",
                                  header=False,
                                  index=False,
                                  )
        
        return super().setUp()

    def tearDown(self) -> None:
        if os.path.exists(self.__test_dir):
            shutil.rmtree(self.__test_dir)

        return super().tearDown()
    
    def get_simple_args(self):
        return argparse.Namespace(reps=["sample_data/ENCFF993VCR.pl.bigWig", "sample_data/ENCFF182TPF.mn.bigWig"],
                                  opath=os.path.join(self.__test_dir, "test_out.bw"),
                                  chrom_size=self.__chrom_size_path,
                                  verbose=False,
                                  )

    def test_main(self):
        args = self.get_simple_args()

        MergeMultiBw.main(args)

        reference_output_path = os.path.join(self.__test_dir, "ref_out.bw")
        merge_bw_files(rep1_bw_path=args.reps[0],
                       rep2_bw_path=args.reps[1],
                       chroms=self.chrom_size_df['chrom'].values,
                       chrom_sizes=self.chrom_size_df['size'].values,
                       opath=reference_output_path, 
                       verbose=False,
                       )
        result_bw = pyBigWig.open(args.opath)

        self.assertFalse(result_bw.intervals("chrFake"))

        ref_bw = pyBigWig.open(reference_output_path)

        for chrom in ["chr3", "chr4"]:
            self.assertEqual(result_bw.intervals(chrom), ref_bw.intervals(chrom))

        result_bw.close()
        ref_bw.close()

