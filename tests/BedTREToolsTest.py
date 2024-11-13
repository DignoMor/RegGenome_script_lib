
import argparse
import unittest
import shutil
import sys
import os

sys.path.append("scripts")
from scripts.bedtretools import BedTRETools
from RGTools.BedTable import BedTable6

class BedTREToolsTest(unittest.TestCase):
    def setUp(self):
        self.__test_dir = "BedTREToolsTest_temp"

        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        self.__sample_bedtre_path = "sample_data/ENCFF819CPA.K562_PROcap.bidirectionalPeak.bedtre"

        return super().setUp()
    
    def tearDown(self):
        shutil.rmtree(self.__test_dir)
        return super().tearDown()
    
    def get_default_args(self):
        args = argparse.Namespace(inpath=self.__sample_bedtre_path,
                                  opath=os.path.join(self.__test_dir, "test_output.bedtre"),
                                  tss_padding=60, 
                                  )
        return args
    
    def test_get_tss_defined_boundry(self):
        args = self.get_default_args()
        
        BedTRETools.main_get_tss_defined_boundry(args)

        output_bt = BedTable6()

        output_bt.load_from_file(args.opath)

        self.assertEqual(output_bt.get_region_names()[0], 
                         "chr1_634027_634244", 
                         )

        self.assertEqual(output_bt.get_chrom_names()[0], 
                         "chr1", 
                         )
        
        self.assertEqual(output_bt.get_start_locs()[0], 
                         633994, 
                         )

        self.assertEqual(output_bt.get_end_locs()[0], 
                         634303, 
                         )
