
import os

import sys
import shutil
import argparse
import pyBigWig
import unittest

import pandas as pd
import numpy as np

sys.path.append("scripts")
from scripts.sample_bw import SampleBw

class SampleBwTest(unittest.TestCase):
    def setUp(self):
        self.__test_dir = "SampleBwTest_temp_data"

        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        self.__chrom_size_path = os.path.join(self.__test_dir, "chrom.size")

        chrom_size_df = pd.DataFrame({"chrom": ["chr1", "chr2", "chr3", "chr4", "chrFake"],
                                      "size": [248956422, 242193529, 198295559, 190214555, 1000000],
                                      },
                                     columns=["chrom", "size"],
                                     )

        chrom_size_df.to_csv(self.__chrom_size_path, 
                             sep="\t", 
                             header=False, 
                             index=False, 
                             )

        return super().setUp()

    def tearDown(self):
        if os.path.exists(self.__test_dir):
            shutil.rmtree(self.__test_dir)

        return super().tearDown()

    def get_simple_args(self):
        return argparse.Namespace(inpath="sample_data/ENCFF993VCR.pl.bigWig",
                                  sample_rate=0.6,
                                  chrom_size=self.__chrom_size_path,
                                  seed=123,
                                  opath=os.path.join(self.__test_dir, "test_out.bw"),
                                  )

    def test_main(self):
        args = self.get_simple_args()
        SampleBw.main(args)

        input_bw = pyBigWig.open(args.inpath)
        output_bw = pyBigWig.open(args.opath)

        self.assertTrue(len(output_bw.intervals("chr1")) < len(input_bw.intervals("chr1")))
        sampled_ratio = np.sum([e[2] for e in output_bw.intervals("chr2")]) / np.sum([e[2] for e in input_bw.intervals("chr2")])
        self.assertAlmostEqual(sampled_ratio, args.sample_rate, places=1)
        self.assertEqual(output_bw.intervals("chr1")[1000][2], 1)
        output_bw.close()

        args.sample_rate = 0.9
        args.opath = os.path.join(self.__test_dir, "test_out_2.bw")
        SampleBw.main(args)

        output_bw = pyBigWig.open(args.opath)

        self.assertTrue(len(output_bw.intervals("chr1")) < len(input_bw.intervals("chr1")))
        sampled_ratio = np.sum([e[2] for e in output_bw.intervals("chr2")]) / np.sum([e[2] for e in input_bw.intervals("chr2")])
        self.assertAlmostEqual(sampled_ratio, args.sample_rate, places=1)
        self.assertEqual(output_bw.intervals("chr1")[1000][2], 4)
        output_bw.close()

        input_bw.close()
