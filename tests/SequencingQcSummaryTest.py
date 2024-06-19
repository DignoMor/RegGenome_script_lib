
import os
import sys
import shutil
import argparse
import unittest

import pandas as pd

sys.path.append('scripts')

from scripts.sequencing_qc_summary import SequencingQcSummary

class SequencingQcSummaryTest(unittest.TestCase):
    def setUp(self) -> None:
        self.__temp_dir = "temp_sequencing_qc_summary_test"
        self.__sample_data_dir = "sample_data"

        if not os.path.exists(self.__temp_dir):
            os.makedirs(self.__temp_dir)

        self.__flagstat_file_list = [os.path.join(self.__sample_data_dir, "sample.bam_flagstat.txt")]
        self.__fastp_json_file_list = [os.path.join(self.__sample_data_dir, "fastp.json")]
        self.__sample_list = ["sample1"]
        return super().setUp()
    
    def tearDown(self) -> None:

        if os.path.exists(self.__temp_dir):
            shutil.rmtree(self.__temp_dir)

        return super().tearDown()
    
    def test_read_flagstat(self):
        flagstat_path = self.__flagstat_file_list[0]
        info_dict = SequencingQcSummary.read_flagstat(flagstat_path)

        self.assertEqual(info_dict["total_reads"], 64495803)
        self.assertEqual(info_dict["mapped_reads"], 64495803)
        self.assertEqual(info_dict["properly_paired"], 45537762)
    
    def test_read_fastp_json(self):
        fastp_json_path = self.__fastp_json_file_list[0]
        info_dict = SequencingQcSummary.read_fastp_json(fastp_json_path)

        self.assertEqual(info_dict["total_reads"], 16763944)
        self.assertEqual(info_dict["passed_filter_reads"], 16034314)
        self.assertEqual(info_dict["low_quality_reads"], 695022)
        self.assertEqual(info_dict["too_many_N_reads"], 14740)
        self.assertEqual(info_dict["too_short_reads"], 19868)
        self.assertEqual(info_dict["too_long_reads"], 0)
        self.assertEqual(info_dict["percentage_passed_filter"], 16034314 / 16763944)
        self.assertEqual(info_dict["percentage_low_quality"], 695022 / 16763944)
        self.assertEqual(info_dict["percentage_too_many_N"], 14740 / 16763944)
        self.assertEqual(info_dict["percentage_too_short"], 19868 / 16763944)
        self.assertEqual(info_dict["percentage_too_long"], 0.0)

    def test_main(self):
        args = argparse.Namespace(
            sample = self.__sample_list,
            flagstat_path = self.__flagstat_file_list,
            fastp_json_path = self.__fastp_json_file_list,
            opath =  os.path.join(self.__temp_dir, "sequencing_qc_summary.txt"), 
            flagstat_fields = ["total_reads", "mapped_reads"], 
            fastp_fields = ["total_reads", "passed_filter_reads", "percentage_passed_filter"], 
            )
        
        SequencingQcSummary.main(args)
        output_df = pd.read_csv(args.opath, 
                                sep="\t", 
                                )
        
        self.assertEqual(output_df.shape[0], 1)
        self.assertEqual(output_df.shape[1], 6)
        self.assertEqual(output_df.loc[0, "sample"], "sample1")
        self.assertEqual(output_df.loc[0, "flagstat_total_reads"], 64495803)
        self.assertEqual(output_df.loc[0, "flagstat_mapped_reads"], 64495803)
        self.assertEqual(output_df.loc[0, "fastp_total_reads"], 16763944)
        self.assertEqual(output_df.loc[0, "fastp_passed_filter_reads"], 16034314)
        self.assertEqual(output_df.loc[0, "fastp_percentage_passed_filter"], 16034314 / 16763944)

        SequencingQcSummary.main(args)

if __name__ == "__main__":
    unittest.main()

