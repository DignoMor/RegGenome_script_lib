
import unittest
import argparse
import shutil
import sys
import os

import pandas as pd

from scripts.RGTools.BedTable import BedTable3, BedTable6, BedTable6Plus

sys.path.append("scripts")
from scripts.count_bw_sig import main

class CountBwSigTest(unittest.TestCase):
    def setUp(self) -> None:
        self.__bw_pls = ["sample_data/ENCFF993VCR.pl.bigWig"]
        self.__bw_mns = ["sample_data/ENCFF182TPF.mn.bigWig"]
        self.__temp_dir = "CountBwSigTest_temp_data"

        if not os.path.exists(self.__temp_dir):
            os.makedirs(self.__temp_dir)

        self.__bed3_path = os.path.join(self.__temp_dir, "test.bed3")
        self.__bed6_path = os.path.join(self.__temp_dir, "test.bed6")
        self.__bed6gene_path = os.path.join(self.__temp_dir, "test.bed6gene")

        region_df = pd.DataFrame({"chr": ["chr6", "chr14", "chr17"],
                                  "start_loc": [170553801, 75278325, 45894026],
                                  "end_loc": [170554802, 75279326, 45895027],
                                  "name": ["TBP", "FOS", "MAPT"],
                                  "score": [1, 2, 3],
                                  "strand": ["+", "-", "+"],
                                  "gene_name": ["gene1", "gene2", "gene3"],
                                  }, 
                                  )

        # Write regions in different bed format
        bedtable3 = BedTable3()
        bedtable3.load_from_dataframe(region_df[["chr", "start_loc", "end_loc"]], 
                                      column_map={"chrom": "chr", 
                                                  "start": "start_loc", 
                                                  "end": "end_loc", 
                                                  },
                                      )
        bedtable3.write(self.__bed3_path)

        bedtable6 = BedTable6()
        bedtable6.load_from_dataframe(region_df[["chr", "start_loc", "end_loc", "name", "score", "strand"]], 
                                      column_map={"chrom": "chr", 
                                                  "start": "start_loc", 
                                                  "end": "end_loc", 
                                                  "name": "name", 
                                                  "score": "score", 
                                                  "strand": "strand", 
                                                  },
                                      )
        bedtable6.write(self.__bed6_path)

        bedtable6gene = BedTable6Plus(extra_column_names=["gene_name"],
                                      extra_column_dtype=[str],
                                      )
        bedtable6gene.load_from_dataframe(region_df, 
                                          column_map={"chrom": "chr", 
                                                      "start": "start_loc", 
                                                      "end": "end_loc", 
                                                      "name": "name", 
                                                      "score": "score", 
                                                      "strand": "strand", 
                                                      "gene_name": "gene_name", 
                                                      },
                                          )
        bedtable6gene.write(self.__bed6gene_path)

        #TODO: Add TRE bed tests

        return super().setUp()

    def tearDown(self) -> None:
        if os.path.exists(self.__temp_dir):
            shutil.rmtree(self.__temp_dir)

        return super().tearDown()

    def test_main_bed3_input(self):
        job_name = "test_bed3_input"
        args = argparse.Namespace(job_name=job_name,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed3",
                                  region_file_path=self.__bed3_path,
                                  single_bw=False,
                                  ignore_strandness=False,
                                  opath=self.__temp_dir,
                                  l_pad = 0, 
                                  r_pad = 0, 
                                  min_len_after_padding=1,
                                  method_resolving_invalid_padding="raise", 
                                  output_type="raw_count",
                                  )
        main(args)

        # Test count output
        count_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".count.csv"), 
                               index_col=0,
                               )
        
        self.assertEqual(count_df.shape, (3, 1))
        self.assertEqual(count_df.loc["chr6_170553801_170554802", "sample1"], 379)

        # test region info output

        region_info_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".region_info.csv"), 
                                     index_col=0,
                                     )
        
        self.assertEqual(region_info_df.shape, (3, 6))
        self.assertEqual(region_info_df.loc["chr6_170553801_170554802", "name"], ".")
    
    def test_main_bed3_input_single_bw(self):
        job_name = "test_bed3_input_single_bw"
        args = argparse.Namespace(job_name=job_name,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed3",
                                  region_file_path=self.__bed3_path,
                                  single_bw=True,
                                  ignore_strandness=True,
                                  opath=self.__temp_dir,
                                  l_pad = 0, 
                                  r_pad = 0, 
                                  min_len_after_padding=1,
                                  method_resolving_invalid_padding="raise", 
                                  output_type="raw_count",
                                  )
        main(args)

        # Test count output
        count_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".count.csv"), 
                               index_col=0,
                               )
        
        self.assertEqual(count_df.shape, (3, 1))
        self.assertEqual(count_df.loc["chr6_170553801_170554802", "sample1"], 348)

        # test region info output

        region_info_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".region_info.csv"), 
                                     index_col=0,
                                     )
        
        self.assertEqual(region_info_df.shape, (3, 6))
        self.assertEqual(region_info_df.loc["chr6_170553801_170554802", "name"], ".")

    def test_main_bed6_input(self):
        job_name = "test_bed6_input"
        args = argparse.Namespace(job_name=job_name,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed6",
                                  region_file_path=self.__bed6_path,
                                  single_bw=False,
                                  ignore_strandness=False,
                                  opath=self.__temp_dir,
                                  l_pad = 0, 
                                  r_pad = 0, 
                                  min_len_after_padding=1,
                                  method_resolving_invalid_padding="raise", 
                                  output_type="raw_count",
                                  )
        
        main(args)

        # Test count output
        count_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".count.csv"), 
                               index_col=0,
                               )
        
        self.assertEqual(count_df.shape, (3, 1))

        self.assertEqual(count_df.loc["chr6_170553801_170554802", "sample1"], 348)
    
    def test_main_bed6gene_input(self):
        job_name = "test_bed6gene_input"
        args = argparse.Namespace(job_name=job_name,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed6gene",
                                  region_file_path=self.__bed6gene_path,
                                  single_bw=False,
                                  ignore_strandness=False,
                                  opath=self.__temp_dir,
                                  l_pad = 0, 
                                  r_pad = 0, 
                                  min_len_after_padding=1,
                                  method_resolving_invalid_padding="raise", 
                                  output_type="raw_count",
                                  )
        main(args)
        # Test count output
        count_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".count.csv"), 
                               index_col=0,
                               )
        
        self.assertEqual(count_df.shape, (3, 1))

        self.assertEqual(count_df.loc["chr6_170553801_170554802", "sample1"], 348)

        # Test region info output
        region_info_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".region_info.csv"),
                                     index_col=0,
                                     )
        self.assertEqual(region_info_df.shape, (3, 7))
        self.assertEqual(region_info_df.loc["chr6_170553801_170554802", "gene_name"], "gene1")

    def test_main_strandness_handle(self):
        # Test ignore strandness for bed6
        job_name_bed6 = "test_bed6_ignore_strandness"

        args = argparse.Namespace(job_name=job_name_bed6,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed6",
                                  region_file_path=self.__bed6_path,
                                  single_bw=False,
                                  ignore_strandness=True,
                                  opath=self.__temp_dir,
                                  l_pad = 0, 
                                  r_pad = 0, 
                                  min_len_after_padding=1,
                                  method_resolving_invalid_padding="raise", 
                                  output_type="raw_count",
                                  )
        
        main(args)

        job_name_bed3 = "test_bed3"

        args = argparse.Namespace(job_name=job_name_bed3,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed3",
                                  region_file_path=self.__bed3_path,
                                  single_bw=False,
                                  ignore_strandness=True,
                                  opath=self.__temp_dir,
                                  l_pad = 0, 
                                  r_pad = 0, 
                                  min_len_after_padding=1,
                                  method_resolving_invalid_padding="raise", 
                                  output_type="raw_count",
                                  )

        main(args)

        bed6_result = pd.read_csv(os.path.join(self.__temp_dir, job_name_bed6 + ".count.csv"), 
                                  index_col=0,
                                  )
                                  
        bed3_result = pd.read_csv(os.path.join(self.__temp_dir, job_name_bed3 + ".count.csv"), 
                                  index_col=0,
                                  )

        self.assertTrue((bed6_result["sample1"] == bed3_result["sample1"]).all())

    def test_padding(self):
        job_name = "test_padding"
        args = argparse.Namespace(job_name=job_name,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed6",
                                  region_file_path=self.__bed6_path,
                                  single_bw=False,
                                  ignore_strandness=False,
                                  opath=self.__temp_dir,
                                  l_pad = -100, 
                                  r_pad = -100, 
                                  min_len_after_padding=1,
                                  method_resolving_invalid_padding="raise", 
                                  output_type="raw_count",
                                  )
        
        main(args)

        # Test count output
        count_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".count.csv"), 
                               index_col=0,
                               )
        
        self.assertEqual(count_df.shape, (3, 1))

        self.assertEqual(count_df.loc["chr6_170553801_170554802", "sample1"], 344)

    def test_resolving_invalid_padding(self):
        job_name = "test_resolving_invalid_padding"
        args = argparse.Namespace(job_name=job_name,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed6",
                                  region_file_path=self.__bed6_path,
                                  single_bw=False,
                                  ignore_strandness=False,
                                  opath=self.__temp_dir,
                                  l_pad = 0, 
                                  r_pad = 0, 
                                  min_len_after_padding=10000,
                                  method_resolving_invalid_padding="raise", 
                                  output_type="raw_count",
                                  )
        
        with self.assertRaises(Exception):
            main(args)
        
        args.min_len_after_padding = 1
        args.l_pad= -1000
        args.r_pad= -1000
        with self.assertRaises(Exception):
            main(args)

        args.method_resolving_invalid_padding = "fallback"
        main(args)

        # Test count output
        count_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".count.csv"), 
                               index_col=0,
                               )
        
        self.assertEqual(count_df.shape, (3, 1))

        self.assertEqual(count_df.loc["chr6_170553801_170554802", "sample1"], 348)
    
    def test_output_type(self):
        job_name = "test_output_type"
        args = argparse.Namespace(job_name=job_name,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed6",
                                  region_file_path=self.__bed6_path,
                                  single_bw=False,
                                  ignore_strandness=False,
                                  opath=self.__temp_dir,
                                  l_pad = -250, 
                                  r_pad = -250,
                                  min_len_after_padding=1,
                                  method_resolving_invalid_padding="raise", 
                                  output_type="RPK",
                                  )
        
        main(args)

        # Test count output
        count_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".RPK.csv"), 
                               index_col=0,
                               )
        
        self.assertEqual(count_df.shape, (3, 1))

        self.assertTrue(count_df.loc["chr6_170553801_170554802", "sample1"] - 680.64 < 0.01)
