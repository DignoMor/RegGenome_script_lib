#!/usr/bin/env python3

import argparse
import unittest
import shutil
import os

import pandas as pd

from scripts.RGTools.BedTable import BedTable3, BedTable6
from scripts.subset_bed import main

class SubsetBedTest(unittest.TestCase):
    def setUp(self) -> None:
        self.__data_dir = "SubsetBedTest_temp_data"
        self.__sample_data_dir = "tests/sample_data"
        self.__sample_gwas_bed6 = os.path.join(self.__sample_data_dir, "sample.gwas_summary.bed6")
        self.__sample_annot_bed3 = os.path.join(self.__sample_data_dir, "sample.annot.bed3")

        if not os.path.exists(self.__data_dir):
            os.makedirs(self.__data_dir)

        return super().setUp()

    def tearDown(self) -> None:
        if os.path.exists(self.__data_dir):
            shutil.rmtree(self.__data_dir)
        return super().tearDown()

    def test_subset_bed(self):
        # test for annot bed3 files
        args = argparse.Namespace(input_bed_path=self.__sample_annot_bed3,
                                  output_bed_path=os.path.join(self.__data_dir, "subset.bed3"),
                                  chr_name="chr2",
                                  start_loc=127100000,
                                  end_loc=127120000,
                                  file_type="bed3", 
                                  )

        main(args)
        
        output_bed_table3 = BedTable3()
        output_bed_table3.load_from_file(args.output_bed_path)

        ans_table = pd.DataFrame({"chr":["chr2", "chr2"], "start":[127105800, 127106884], "end":[127106242, 127107434]})

        self.assertArrayEqual(output_bed_table3.to_dataframe().values, 
                              ans_table.values, 
                              )

        # test for gwas bed6 files

        args = argparse.Namespace(input_bed_path=self.__sample_gwas_bed6,
                                  output_bed_path=os.path.join(self.__data_dir, "subset.bed6"),
                                  chr_name="chr2",
                                  start_loc=127100000,
                                  end_loc=127101000,
                                  file_type="bed6", 
                                  )
        
        main(args)

        output_bed_table6 = BedTable6()
        output_bed_table6.load_from_file(args.output_bed_path)

        snp_locs = [127100135, 127100370, 127100560, 127100633, 127100891, 127100953]

        ans_table = pd.DataFrame({"chr":["chr2", "chr2", "chr2", "chr2", "chr2", "chr2"], 
                                  "start":snp_locs, 
                                  "end":[loc + 1 for loc in snp_locs], 
                                  "name":["rs188131231", "rs734115", "rs12618318", "rs185842969", "rs190671778", "rs550084979"], 
                                  "score":[0.8315, 0.007222, 0.2419, 0.9132, 0.006013, 0.5691], 
                                  "strand":["+", "+", "+", "+", "+", "+"], 
                                  }, 
                                  )
        
        self.assertArrayEqual(output_bed_table6.to_dataframe().values,
                              ans_table.values,
                              )

    def assertArrayEqual(self, arr1, arr2):

        self.assertEqual(len(arr1), len(arr2))
        self.assertTrue((arr1 == arr2).all())

if __name__ == "__main__":
    unittest.main()
