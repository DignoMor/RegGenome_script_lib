
import unittest
import shutil
import os

import numpy as np
import pandas as pd

from scripts.RGTools.BedTable import BedTable3

class TestBedTable3(unittest.TestCase):
    def setUp(self) -> None:
        self.__data_dir = "bedTableTest_temp_data"
        if not os.path.exists(self.__data_dir):
            os.mkdir(self.__data_dir)
        self.__data_file = os.path.join(self.__data_dir, "test.bed")

        self.__data_df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "start": [1, 8, 3, 4],
            "end": [5, 12, 7, 8],
        })

        self.__data_df.to_csv(self.__data_file, 
                              sep="\t", 
                              header=False, 
                              index=False, 
                              )
        
        return super().setUp()
    
    def tearDown(self) -> None:
        if os.path.exists(self.__data_dir):
            shutil.rmtree(self.__data_dir)

        return super().tearDown()
    
    def test_load_from_file(self):
        bed_table = BedTable3()
        bed_table.load_from_file(self.__data_file)

        self.assertEqual(bed_table.to_dataframe().values, 
                         self.__data_df.values,
                         )
        
    def test_load_from_dataframe(self):
        bed_table = BedTable3()
        bed_table.load_from_dataframe(self.__data_df.copy())

        self.assertEqual(bed_table.to_dataframe().values, 
                         self.__data_df.values,
                         )

    def test_apply_logical_filter(self):
        logical_array = [True, False, True, False]
        logical_array = np.array(logical_array)

        bed_table = self.__init_test_bed_table()

        filtered_bed_table = bed_table.apply_logical_filter(logical_array)

        self.assertEqual(filtered_bed_table.to_dataframe().values, 
                         self.__data_df.loc[logical_array].values, 
                         )

    def test_region_subset(self):
        bed_table = self.__init_test_bed_table()

        chrom = "chr2"
        lboundary = 2
        rboundary = 7

        subset_bed_table = bed_table.region_subset(chrom=chrom,
                                                   start=lboundary, 
                                                   end=rboundary,
                                                   )

        self.assertEqual(subset_bed_table.to_dataframe().values,
                         self.__data_df.loc[(self.__data_df["chrom"] == chrom) & \
                                            (self.__data_df["start"] >= lboundary) & \
                                             self.__data_df["end"] <= rboundary].values,
                         )

    def test_to_dataframe(self):
        bed_table = self.__init_test_bed_table()

        exported_dataframe = bed_table.to_dataframe()

        self.assertEqual(exported_dataframe.values,
                         self.__data_df.values,
                         )
        
        # Test immutability
        exported_dataframe.loc[0, 0] = "chr3"
        self.assertTrue(self.__mutate_test_helper(bed_table))

    def test_write(self):
        bed_table = self.__init_test_bed_table()

        write_path = os.path.join(self.__data_dir, "write_test.bed")
        bed_table.write(write_path)

        new_bed_table = BedTable3()
        new_bed_table.load_from_file(write_path)

        self.assertEqual(new_bed_table.to_dataframe().values,
                         bed_table.to_dataframe().values,
                         )

    def test_get_chrom_names(self):
        bed_table = self.__init_test_bed_table()

        chrom_names = bed_table.get_chrom_names()

        self.assertEqual(chrom_names,
                         self.__data_df["chrom"].values,
                         )

    def test_get_start_locs(self):
        bed_table = self.__init_test_bed_table()

        start_locs = bed_table.get_start_locs()

        self.assertEqual(start_locs,
                         self.__data_df["start"].values,
                         )

    def test_get_end_locs(self):
        bed_table = self.__init_test_bed_table()

        end_locs = bed_table.get_end_locs()

        self.assertEqual(end_locs, 
                         self.__data_df["end"].values,
                         )

    def test_iter_regions(self):
        bed_table = self.__init_test_bed_table()
        start_locs = np.array([r["start"] for r in bed_table.iter_regions()])

        self.assertEqual(start_locs,
                         bed_table.get_start_locs(),
                         )

    def __mutate_test_helper(self, bedtable: BedTable3) -> bool:
        '''
        Helper function to test the immutability of the class.
        Return True if the passed object is not mutated.
        '''
        return self.__data_df.values == bedtable.to_dataframe().values
    
    def __init_test_bed_table(self) -> BedTable3:
        bed_table = BedTable3()
        bed_table.load_from_file(self.__data_file)

        return bed_table

if __name__ == "__main__":
    unittest.main()
