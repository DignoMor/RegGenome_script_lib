
import unittest
import argparse
import shutil
import sys
import os

import pandas as pd

sys.path.append("scripts")

from scripts.padding_bed import PaddingBed
from RGTools.BedTable import BedTable6, BedTable3

class PaddingBedTest(unittest.TestCase):
    def setUp(self):
        self.__test_dir = "PaddingBedTest_temp"
        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        self.__test_bed6_path = os.path.join(self.__test_dir, "test_input.bed6")
        self.__test_bed3_path = os.path.join(self.__test_dir, "test_input.bed3")
        self.__test_opath = os.path.join(self.__test_dir, "test_output.bed")

        input_bt3 = BedTable3()
        input_bt6 = BedTable6()
        input_df = pd.DataFrame({"chrom": ["chr1", "chr2", "chr3"],
                                 "start": [10001, 11001, 12001],
                                 "end": [10002, 11002, 12002],
                                 "name": ["gene1", "gene2", "gene3"],
                                 "score": [0, 0, 0],
                                 "strand": ["+", "-", "-"]
                                 }, 
                                columns=["chrom", "start", "end", "name", "score", "strand"], 
                                )

        input_bt6.load_from_dataframe(input_df)
        input_bt3.load_from_dataframe(input_df[["chrom", "start", "end"]])
    
        input_bt6.write(self.__test_bed6_path)
        input_bt3.write(self.__test_bed3_path)

        return super().setUp()
    
    def tearDown(self):
        shutil.rmtree(self.__test_dir)
        return super().tearDown()

    def get_simple_params(self):
        args = argparse.Namespace(inpath=self.__test_bed6_path,
                                  opath=self.__test_opath, 
                                  upstream_pad=400,
                                  downstream_pad=500, 
                                  input_file_type="bed6", 
                                  ignore_strand=False,
                                  )
        return args
    
    def test_simple_padding(self):
        args = self.get_simple_params()

        PaddingBed.main(args)
        output_bt = BedTable6()
        output_bt.load_from_file(self.__test_opath)
        output_df = output_bt.to_dataframe()

        self.assertEqual(output_df.shape[0], 3)
        self.assertEqual(output_df.loc[0, "start"], 9601)
        self.assertEqual(output_df.loc[0, "end"], 10502)
        self.assertEqual(output_df.loc[1, "start"], 10501)
        self.assertEqual(output_df.loc[1, "end"], 11402)

    def test_bed3_input(self):
        args = self.get_simple_params()
        args.input_file_type = "bed3"
        args.inpath = self.__test_bed3_path
        args.ignore_strand = True

        PaddingBed.main(args)

        output_bt = BedTable3()
        output_bt.load_from_file(self.__test_opath)
        output_df = output_bt.to_dataframe()

        self.assertEqual(output_df.shape[0], 3)
        self.assertEqual(output_df.loc[0, "start"], 9601)
        self.assertEqual(output_df.loc[0, "end"], 10502)
        self.assertEqual(output_df.loc[1, "start"], 10601)
        self.assertEqual(output_df.loc[1, "end"], 11502)


