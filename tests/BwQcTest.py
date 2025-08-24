
import unittest
import argparse
import shutil
import sys
import os

import pandas as pd

sys.path.append("scripts")

from scripts.bw_qc import BwQc


class BwQcTest(unittest.TestCase):
    def setUp(self) -> None:
        self.__temp_dir = "BwQcTest_temp_data"

        self.__bw_pl = "sample_data/ENCFF993VCR.pl.bigWig"
        self.__bw_mn = "sample_data/ENCFF182TPF.mn.bigWig"
        self.__chrom_sizes = os.path.join(self.__temp_dir, "hg38.chrom.sizes")

        if not os.path.exists(self.__temp_dir):
            os.makedirs(self.__temp_dir)

        pd.DataFrame({"chr": ["chr1", "chr2", "chr3", "chr4", "chr5"],
                      "size": [248956422, 242193529, 198295559, 190214555, 181538259],
                      }).to_csv(self.__chrom_sizes, sep="\t", index=False, header=False)
    
    def tearDown(self) -> None:
        if os.path.exists(self.__temp_dir):
            shutil.rmtree(self.__temp_dir)

    def test_calculate_total_counts(self):
        args = argparse.Namespace(bw_pl=self.__bw_pl,
                                  bw_mn=self.__bw_mn,
                                  chrom_sizes=self.__chrom_sizes,
                                  statistic="total_counts",
                                  opath=os.path.join(self.__temp_dir, "total_counts.txt"),
                                  )
        BwQc.main(args)

        with open(os.path.join(self.__temp_dir, "total_counts.txt"), "r") as f:
            self.assertEqual(f.read(), "3957852.00")

        args.bw_mn = None
        BwQc.main(args)

        with open(os.path.join(self.__temp_dir, "total_counts.txt"), "r") as f:
            self.assertEqual(f.read(), "1950600.00")

    