
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
        self.__bed6plus_path = os.path.join(self.__temp_dir, "test.bed6plus")

        region_df = pd.DataFrame({"chr": ["chr6", "chr14", "chr17"],
                                  "start_loc": [170553801, 75278325, 45894026],
                                  "end_loc": [170554802, 75279326, 45895027],
                                  "name": ["TBP", "FOS", "MAPT"],
                                  "score": [1, 2, 3],
                                  "strand": ["+", "-", "+"],
                                  "extra": ["extra1", "extra2", "extra3"],
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

        bedtable6plus = BedTable6Plus(extra_column_names=["extra"],
                                      extra_column_dtype=[str],
                                      )
        bedtable6plus.load_from_dataframe(region_df, 
                                          column_map={"chrom": "chr", 
                                                      "start": "start_loc", 
                                                      "end": "end_loc", 
                                                      "name": "name", 
                                                      "score": "score", 
                                                      "strand": "strand", 
                                                      "extra": "extra", 
                                                      },
                                          )
        bedtable6plus.write(self.__bed6plus_path)

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
                                  ignore_strandness=False,
                                  opath=self.__temp_dir,
                                  region_padding="0,0", 
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

    def test_main_bed6_input(self):
        job_name = "test_bed6_input"
        args = argparse.Namespace(job_name=job_name,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed6",
                                  region_file_path=self.__bed6_path,
                                  ignore_strandness=False,
                                  opath=self.__temp_dir,
                                  region_padding="0,0", 
                                  output_type="raw_count",
                                  )
        
        main(args)

        # Test count output
        count_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".count.csv"), 
                               index_col=0,
                               )
        
        self.assertEqual(count_df.shape, (3, 1))

        self.assertEqual(count_df.loc["chr6_170553801_170554802", "sample1"], 348)
    
    def test_main_strandness_handle(self):
        # Test ignore strandness for bed6
        job_name_bed6 = "test_bed6_ignore_strandness"

        args = argparse.Namespace(job_name=job_name_bed6,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed6",
                                  region_file_path=self.__bed6_path,
                                  ignore_strandness=True,
                                  opath=self.__temp_dir,
                                  region_padding="0,0", 
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
                                  ignore_strandness=True,
                                  opath=self.__temp_dir,
                                  region_padding="0,0", 
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
                                  ignore_strandness=False,
                                  opath=self.__temp_dir,
                                  region_padding="-100,-100", 
                                  output_type="raw_count",
                                  )
        
        main(args)

        # Test count output
        count_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".count.csv"), 
                               index_col=0,
                               )
        
        self.assertEqual(count_df.shape, (3, 1))

        self.assertEqual(count_df.loc["chr6_170553801_170554802", "sample1"], 344)
    
    def test_output_type(self):
        job_name = "test_output_type"
        args = argparse.Namespace(job_name=job_name,
                                  sample_names=["sample1"],
                                  bw_pls=self.__bw_pls, 
                                  bw_mns=self.__bw_mns,
                                  region_file_type="bed6",
                                  region_file_path=self.__bed6_path,
                                  ignore_strandness=False,
                                  opath=self.__temp_dir,
                                  region_padding="-250,-250", 
                                  output_type="RPK",
                                  )
        
        main(args)

        # Test count output
        count_df = pd.read_csv(os.path.join(self.__temp_dir, job_name + ".RPK.csv"), 
                               index_col=0,
                               )
        
        self.assertEqual(count_df.shape, (3, 1))

        self.assertTrue(count_df.loc["chr6_170553801_170554802", "sample1"] - 680.64 < 0.01)
