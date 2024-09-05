
import os
import sys
import shutil
import argparse
import unittest
import pyBigWig

import pandas as pd

sys.path.append("scripts")
from scripts.filter_bw import FilterBw

class FilterBwTest(unittest.TestCase):
    def setUp(self) -> None:
        self.__test_dir = "filter_bw_test"

        self.__chrom_size_path = os.path.join(self.__test_dir, "chrom.size")

        chrom_size_df = pd.DataFrame({"chrom": ["chr1", "chr2", "chr3", "chr4", "chrFake"],
                                      "size": [248956422, 242193529, 198295559, 190214555, 1000000],
                                      },
                                     columns=["chrom", "size"],
                                     )
        
        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        chrom_size_df.to_csv(self.__chrom_size_path, 
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
        return argparse.Namespace(inpath="sample_data/ENCFF993VCR.pl.bigWig", 
                                  opath="test_out.bw", 
                                  chrom_size=self.__chrom_size_path, 
                                  )
    
    def test_main(self):
        args = self.get_simple_args()

        FilterBw.main(args)

        input_bw = pyBigWig.open(args.inpath)
        result_bw = pyBigWig.open(args.opath)

        self.assertFalse(result_bw.intervals("chrFake"))

        for chrom in ["chr1", "chr2", "chr3", "chr4"]:
            self.assertEqual(result_bw.intervals(chrom), 
                             input_bw.intervals(chrom),
                             )
