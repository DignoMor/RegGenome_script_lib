
import argparse
import unittest
import shutil
import sys
import os

import pandas as pd

sys.path.append("scripts")
from scripts.bed2tssbed import Bed2TSSBED
from scripts.RGTools.BedTable import BedTable3, BedTable6

class Bed2tssbedTest(unittest.TestCase):
    def setUp(self) -> None:
        self.__data_dir = "Bed2tssbedTest_temp_data"

        if not os.path.exists(self.__data_dir):
            os.makedirs(self.__data_dir)

        self.__test_data = pd.DataFrame({"chrom": ["chr1", "chr2", "chr3"],
                                         "start": [100, 200, 300],
                                         "end": [200, 300, 400],
                                         "name": ["gene1", "gene2", "gene3"],
                                         "score": [1, 2, 3],
                                         "strand": ["+", "-", "+"],
                                         })
        
        self.__test_ans = pd.DataFrame({"chrom": ["chr1", "chr2", "chr3"],
                                        "start": [100, 299, 300],
                                        "end": [101, 300, 301],
                                        "name": ["gene1", "gene2", "gene3"],
                                        "score": [1, 2, 3],
                                        "strand": ["+", "-", "+"],
                                        })
        
        self.__test_bed_table = BedTable6()
        self.__test_bed_table.load_from_dataframe(self.__test_data)

        self.__ans_bed_table = BedTable6()
        self.__ans_bed_table.load_from_dataframe(self.__test_ans)

        self.__bed_in = os.path.join(self.__data_dir, "test.bed")
        self.__bed_out = os.path.join(self.__data_dir, "test.tss.bed")
        self.__bed_in_bed6gene = os.path.join(self.__data_dir, "test.bed6gene")
        self.__bed_out_bed6gene = os.path.join(self.__data_dir, "test.tss.bed6gene")
        
        self.__test_bed_table.write(self.__bed_in)

        # test for bed6gene
        self.__test_data_bed6gene = self.__test_data.copy()
        self.__test_ans_bed6gene = self.__test_ans.copy()

        self.__test_data_bed6gene["gene_symbol"] = ["gene1", "gene2", "gene3"]
        self.__test_ans_bed6gene["gene_symbol"] = ["gene1", "gene2", "gene3"]

        self.__test_bt_bed6gene = Bed2TSSBED.BedTable6Gene()
        self.__test_bt_bed6gene.load_from_dataframe(self.__test_data_bed6gene)
        self.__ans_bt_bed6gene = Bed2TSSBED.BedTable6Gene()
        self.__ans_bt_bed6gene.load_from_dataframe(self.__test_ans_bed6gene)

        self.__test_bt_bed6gene.write(self.__bed_in_bed6gene)
        
        return super().setUp()

    def tearDown(self) -> None:
        if os.path.exists(self.__data_dir):
            shutil.rmtree(self.__data_dir)
        return super().tearDown()
    
    def get_simple_args(self):
        args = argparse.Namespace(bed_in=self.__bed_in,
                                  bed_out=self.__bed_out,
                                  region_file_type="bed6",
                                  )
        
        return args
    
    def test_bed2tssbed(self):
        args = self.get_simple_args()

        Bed2TSSBED.main(args)

        self.__out_bed_table = BedTable6()
        self.__out_bed_table.load_from_file(self.__bed_out)

        self.assertTrue((self.__out_bed_table.to_dataframe().values == self.__ans_bed_table.to_dataframe().values).all())
    
    def test_bed6gene_io(self):
        args = self.get_simple_args()

        args.region_file_type = "bed6gene"
        args.bed_in = self.__bed_in_bed6gene
        args.bed_out = self.__bed_out_bed6gene

        Bed2TSSBED.main(args)

        self.__out_bt_bed6gene = Bed2TSSBED.BedTable6Gene()
        self.__out_bt_bed6gene.load_from_file(args.bed_out)

        self.assertTrue((self.__out_bt_bed6gene.to_dataframe().values == self.__ans_bt_bed6gene.to_dataframe().values).all())
    
if __name__ == "__main__":
    unittest.main()
