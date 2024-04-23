
import argparse
import unittest
import shutil
import os

import pandas as pd

from scripts.bed2tssbed import bed2tssbed
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
        
        self.__test_bed_table.write(self.__bed_in)

        
        return super().setUp()

    def tearDown(self) -> None:
        if os.path.exists(self.__data_dir):
            shutil.rmtree(self.__data_dir)
        return super().tearDown()
    
    def test_bed2tssbed(self):
        args = argparse.Namespace(bed_in=self.__bed_in,
                                  bed_out=self.__bed_out,
                                  window_size="0-0",
                                  )

        bed2tssbed(args)

        self.__out_bed_table = BedTable6()
        self.__out_bed_table.load_from_file(self.__bed_out)

        self.assertTrue((self.__out_bed_table.to_dataframe().values == self.__ans_bed_table.to_dataframe().values).all())
    
if __name__ == "__main__":
    unittest.main()
