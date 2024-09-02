
import os
import sys
import shutil
import argparse
import unittest

from scripts.RGTools.BedTable import BedTable6

sys.path.append("scripts")
from scripts.identify_pausing_region import IdentifyPausingRegion

class IdentifyPausingRegionTest(unittest.TestCase):
    def setUp(self) -> None:
        self.__tss_annotation_path = "sample_data/gencode.v37.protein_coding_gene_body.tss.chr2.bed6"
        self.__bw_pl_path = "sample_data/ENCFF993VCR.pl.bigWig"
        self.__bw_mn_path = "sample_data/ENCFF182TPF.mn.bigWig"

        self.__test_folder = "identify_pausing_region_test"

        if not os.path.exists(self.__test_folder):
            os.makedirs(self.__test_folder)

        return super().setUp()

    def tearDown(self) -> None:
        shutil.rmtree(self.__test_folder)
        return super().tearDown()
    
    def get_simple_args(self):
        return argparse.Namespace(tss_bed=self.__tss_annotation_path,
                                  l_pad=500,
                                  r_pad=499,
                                  bw_pl=self.__bw_pl_path,
                                  bw_mn=self.__bw_mn_path, 
                                  target_size=250,
                                  search_step_size=10,
                                  opath=os.path.join(self.__test_folder, "pausing_region.bed")
                                  )

    def test_output_size(self):

        args = self.get_simple_args()

        tss_bt = BedTable6()
        tss_bt.load_from_file(self.__tss_annotation_path)

        IdentifyPausingRegion.main(args)
        result_bt = BedTable6()
        result_bt.load_from_file(args.opath)

        self.assertEqual(len(tss_bt), len(result_bt))

    def test_output_rpk(self):
        args = self.get_simple_args()

        IdentifyPausingRegion.main(args)
        result_bt = BedTable6()
        result_bt.load_from_file(args.opath)
        result_df = result_bt.to_dataframe()
        result_df.loc[result_df["name"] == "ENSG00000182551.14", "score"] == 728
