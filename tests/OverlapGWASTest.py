
import argparse
import unittest
import shutil
import sys
import os

sys.path.append("scripts")

from scripts.overlap_GWAS import OverlapGWAS

class OverlapGWASTest(unittest.TestCase):
    def setUp(self):
        self.__test_dir = "OverlapGWASTest_temp_data"
        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        self.__sample_region_bed = os.path.join("sample_data/sample.snp_overlapping_tre.chr2.bed")
        self.__GWAS_file_path = os.path.join("sample_data/sample.chr2.snp.bed")

        self.__GWAS_path_by_tiers = []
        self.__GWAS_path_by_tiers.append(os.path.join(self.__test_dir, "sample.chr2.tier1.snp.bed"))
        self.__GWAS_path_by_tiers.append(os.path.join(self.__test_dir, "sample.chr2.tier2.snp.bed"))
        self.__GWAS_path_by_tiers.append(os.path.join(self.__test_dir, "sample.chr2.tier3.snp.bed"))

        gwas_bt = OverlapGWAS.SNPBed()

        gwas_bt.load_from_file(self.__GWAS_file_path)
        gwas_bt.apply_logical_filter(gwas_bt.get_region_scores() < 5e-15).write(self.__GWAS_path_by_tiers[0])
        gwas_bt.apply_logical_filter(gwas_bt.get_region_scores() < 5e-12).write(self.__GWAS_path_by_tiers[1])
        gwas_bt.apply_logical_filter(gwas_bt.get_region_scores() > 5e-12).write(self.__GWAS_path_by_tiers[2])

        return super().setUp()
    
    def tearDown(self):
        if os.path.exists(self.__test_dir):
            shutil.rmtree(self.__test_dir)
        return super().tearDown()
    
    def get_simple_args(self):
        args = argparse.Namespace(job_name="test",
                                  region_path=self.__sample_region_bed, 
                                  opath=os.path.join(self.__test_dir),
                                  SNP_paths=self.__GWAS_path_by_tiers,
                                  region_file_type="bed6",
                                  )
        return args
    
    def test_main(self):
        args = self.get_simple_args()

        OverlapGWAS.main(args)

        output_bt = OverlapGWAS.SNPBed()
        output_bt.load_from_file(os.path.join(args.opath, "{}.tier{:d}.bed".format(args.job_name, 1)))
        self.assertEqual(output_bt.get_start_locs()[0], 127084012)
        self.assertEqual(output_bt.get_end_locs()[1], 127086450)

        output_bt = OverlapGWAS.SNPBed()
        output_bt.load_from_file(os.path.join(args.opath, "{}.tier{:d}.bed".format(args.job_name, 2)))
        self.assertEqual(len(output_bt), 1)
