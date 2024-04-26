#!/usr/bin/env python3

import argparse
import unittest
import shutil
import sys
import os

import pandas as pd

from scripts.RGTools.BedTable import BedTable6, BedTable6Plus

sys.path.append("scripts")
from scripts.gtf_search import main

class GTFSearchTest(unittest.TestCase):
    def setUp(self) -> None:
        self.__data_dir = "GTFSearchTest_temp_data"
        self.__sample_gtf_path = "sample_data/sample.GENCODE.gtf"

        if not os.path.exists(self.__data_dir):
            os.makedirs(self.__data_dir)

        return super().setUp()
    
    def tearDown(self) -> None:
        if os.path.exists(self.__data_dir):
            shutil.rmtree(self.__data_dir)

        return super().tearDown()
    
    def test_gtf_search(self):

        # test the result of gtf_search
        # output all the protein transcript annotations on chr2
        args = argparse.Namespace(gtf_path=self.__sample_gtf_path,
                                  bed_out=os.path.join(self.__data_dir, "gtf_search.bed"), 
                                  general_feature_key_value_pair=["chr_name==chr2", "feature_type==transcript"],
                                  additional_feature_key_value_pair=["gene_type==protein_coding"],
                                  extra_col_general_feature=[],
                                  extra_col_additional_feature=[],
                                  )

        main(args)
        
        output_bed_table = BedTable6()
        output_bed_table.load_from_file(args.bed_out)

        self.assertEqual(len(output_bed_table), 14)
        self.assertEqual(output_bed_table.get_region_by_index(0)["start"], 127048026)
        self.assertEqual(output_bed_table.get_region_by_index(0)["end"], 127107000)
        self.assertEqual(output_bed_table.get_region_by_index(10)["start"], 127048032)
        self.assertEqual(output_bed_table.get_region_by_index(10)["end"], 127056490)
        self.assertTrue((output_bed_table.get_region_strands() == "-").all())

        # test bed table plus mode
        # Out put all the protein coding transcript annotations on chr2 with gene name as extra column

        args = argparse.Namespace(gtf_path=self.__sample_gtf_path,
                                  bed_out=os.path.join(self.__data_dir, "gtf_search_w_extra_col.bed"), 
                                  general_feature_key_value_pair=["chr_name==chr2", "feature_type==transcript"],
                                  additional_feature_key_value_pair=["gene_type==protein_coding"],
                                  extra_col_general_feature=[],
                                  extra_col_additional_feature=["gene_name"],
                                  )
        main(args)
        
        output_bed_table = BedTable6Plus(extra_column_names=["gene_name"], 
                                         extra_column_dtype=[str], 
                                         )
        output_bed_table.load_from_file(args.bed_out)

        self.assertTrue((output_bed_table.get_region_extra_column("gene_name") == "BIN1").all())

        # test a input and output with no filtering
        args = argparse.Namespace(gtf_path=self.__sample_gtf_path,
                                  bed_out=os.path.join(self.__data_dir, "gtf_search_no_filter.bed"), 
                                  general_feature_key_value_pair=[],
                                  additional_feature_key_value_pair=[],
                                  extra_col_general_feature=[],
                                  extra_col_additional_feature=[],
                                  )
        
        main(args)

        output_bed_table = BedTable6()
        output_bed_table.load_from_file(args.bed_out)

        self.assertEqual(len(output_bed_table), 423)
        