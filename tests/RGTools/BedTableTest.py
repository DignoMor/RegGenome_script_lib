
import unittest
import shutil
import os

import numpy as np
import pandas as pd

from scripts.RGTools.BedTable import BedTable3, BedTable6, BedTable6Plus

class TestBedTable(unittest.TestCase):
    # Public methods for testing BedTable
    def assertArrayEqual(self, arr1, arr2):
        '''
        equal assertion to compare two numpy arrays.
        '''
        if arr1.shape != arr2.shape:
            return self.assertFalse(True, 
                                    msg="Shape mismatch in array comparison: {} != {}".format(arr1.shape, arr2.shape),
                                    )

        return self.assertTrue((arr1 == arr2).all(), 
                               msg="Array mismatch in array comparison: {} != {}".format(arr1, arr2),
                               )

    def mutate_test_helper(self, bedtable: BedTable3) -> bool:
        '''
        Helper function to test the immutability of the class.
        Return True if the passed object is not mutated.
        '''
        return self.assertArrayEqual(self.data_df.values, 
                                     bedtable.to_dataframe().values, 
                                     )
    
    def setUp(self) -> None:
        super().setUp()
        self.data_dir = "bedTableTest_temp_data"
        if not os.path.exists(self.data_dir):
            os.mkdir(self.data_dir)
        self.data_file = os.path.join(self.data_dir, "test.bed")

    def tearDown(self) -> None:
        super().tearDown()
        if os.path.exists(self.data_dir):
            shutil.rmtree(self.data_dir)

class TestBedTable3(TestBedTable):
    def setUp(self) -> None:
        super().setUp()

        self.data_df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "start": [1, 8, 3, 4],
            "end": [5, 12, 7, 8],
        })

        self.data_df.to_csv(self.data_file, 
                              sep="\t", 
                              header=False, 
                              index=False, 
                              )
    
    def tearDown(self) -> None:
        super().tearDown()
    
    def test_load_from_file(self):
        bed_table = BedTable3()
        bed_table.load_from_file(self.data_file)

        self.assertArrayEqual(bed_table.to_dataframe().values, 
                              self.data_df.values,
                              )
        
    def test_load_from_dataframe(self):
        bed_table = BedTable3()
        bed_table.load_from_dataframe(self.data_df.copy())

        self.assertArrayEqual(bed_table.to_dataframe().values, self.data_df.values)

    def test_apply_logical_filter(self):
        logical_array = [True, False, True, False]
        logical_array = np.array(logical_array)

        bed_table = self.__init_test_bed_table()

        filtered_bed_table = bed_table.apply_logical_filter(logical_array)

        self.assertArrayEqual(filtered_bed_table.to_dataframe().values, self.data_df.loc[logical_array].values)

    def test_region_subset(self):
        bed_table = self.__init_test_bed_table()

        chrom = "chr2"
        lboundary = 2
        rboundary = 7

        subset_bed_table = bed_table.region_subset(chrom=chrom,
                                                   start=lboundary, 
                                                   end=rboundary,
                                                   )

        self.assertArrayEqual(subset_bed_table.to_dataframe().values,
                              self.data_df.loc[(self.data_df["chrom"] == chrom) & \
                                               (self.data_df["start"] >= lboundary) & \
                                               (self.data_df["end"] <= rboundary)].values,
                              )

    def test_to_dataframe(self):
        bed_table = self.__init_test_bed_table()

        exported_dataframe = bed_table.to_dataframe()

        self.assertEqual(type(exported_dataframe), type(pd.DataFrame()))
        self.assertArrayEqual(exported_dataframe.values,
                              self.data_df.values,
                              )
        
        # Test immutability
        exported_dataframe.loc[0, 0] = "chr3"
        self.mutate_test_helper(bed_table)

    def test_write(self):
        bed_table = self.__init_test_bed_table()

        write_path = os.path.join(self.data_dir, "write_test.bed")
        bed_table.write(write_path)

        new_bed_table = BedTable3()
        new_bed_table.load_from_file(write_path)

        self.assertArrayEqual(new_bed_table.to_dataframe().values,
                              bed_table.to_dataframe().values,
                              )

    def test_get_chrom_names(self):
        bed_table = self.__init_test_bed_table()

        chrom_names = bed_table.get_chrom_names()

        self.assertArrayEqual(chrom_names, self.data_df["chrom"].values)

    def test_get_start_locs(self):
        bed_table = self.__init_test_bed_table()

        start_locs = bed_table.get_start_locs()

        self.assertArrayEqual(start_locs, self.data_df["start"].values)

    def test_get_end_locs(self):
        bed_table = self.__init_test_bed_table()

        end_locs = bed_table.get_end_locs()

        self.assertArrayEqual(end_locs, self.data_df["end"].values)

    def test_iter_regions(self):
        bed_table = self.__init_test_bed_table()
        start_locs = np.array([r["start"] for r in bed_table.iter_regions()])

        self.assertArrayEqual(start_locs, self.data_df["start"].values)

    def test_sort(self):
        bed_table = BedTable3()
        bed_table.load_from_file(self.data_file)

        bed_table._sort()

        self.assertArrayEqual(bed_table.to_dataframe().values,
                              self.data_df.sort_values(["chrom", "start"]).values,
                              )

    def __init_test_bed_table(self) -> BedTable3:
        bed_table = BedTable3()
        bed_table.load_from_file(self.data_file)

        return bed_table

class TestBedTable6(TestBedTable):
    def setUp(self) -> None:
        super().setUp()
        self.data_df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "start": [1, 8, 3, 4],
            "end": [5, 12, 7, 8],
            "name": ["name1", "name2", "name3", "name4"],
            "score": [0.1, 0.2, 0.3, 0.4],
            "strand": ["+", "-", "+", "-"],
        })

        self.data_df.to_csv(self.data_file, 
                            sep="\t", 
                            header=False, 
                            index=False, 
                            )
        
    def tearDown(self) -> None:
        super().tearDown()
    
    def test_load_from_file(self):
        bed_table = BedTable6()
        bed_table.load_from_file(self.data_file)

        self.assertArrayEqual(bed_table.to_dataframe().values, 
                              self.data_df.values,
                              )
        
    def test_load_from_dataframe(self):
        bed_table = BedTable6()
        bed_table.load_from_dataframe(self.data_df.copy())

        self.assertArrayEqual(bed_table.to_dataframe().values, self.data_df.values)
    
    def test_get_region_names(self):
        bed_table = self.__init_test_bed_table()

        region_names = bed_table.get_region_names()

        self.assertArrayEqual(region_names, self.data_df["name"].values)

    def test_get_region_scores(self):
        bed_table = self.__init_test_bed_table()

        region_scores = bed_table.get_region_scores()

        self.assertArrayEqual(region_scores, self.data_df["score"].values)

    def test_get_region_strands(self):
        bed_table = self.__init_test_bed_table()

        region_strands = bed_table.get_region_strands()

        self.assertArrayEqual(region_strands, self.data_df["strand"].values)

    def __init_test_bed_table(self) -> BedTable6:
        bed_table = BedTable6()
        bed_table.load_from_file(self.data_file)

        return bed_table

class TestBedTable6Plus(TestBedTable):
    def setUp(self) -> None:
        super().setUp()

        self.data_df = pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "start": [1, 8, 3, 4],
            "end": [5, 12, 7, 8],
            "name": ["name1", "name2", "name3", "name4"],
            "score": [0.1, 0.2, 0.3, 0.4],
            "strand": ["+", "-", "+", "-"],
            "extra_str_field": ["extra1", "extra2", "extra3", "extra4"],
            "extra_int_field": [1, 2, 3, 4],
        })

        self.extra_field_names = list(self.data_df.columns[6:])
        self.extra_field_dtype = [str, int]

        self.data_df.to_csv(self.data_file, 
                            sep="\t", 
                            header=False, 
                            index=False, 
                            )

    def tearDown(self) -> None:
        return super().tearDown()

    def test_get_region_extra_column(self):
        bed_table = self.__init_test_bed_table()

        for extra_field_name, extra_field_dtype in zip(self.extra_field_names, self.extra_field_dtype):
            extra_field_data = bed_table.get_region_extra_column(extra_field_name)

            self.assertArrayEqual(extra_field_data, self.data_df[extra_field_name].values)

    def test_load_from_dataframe(self):
        bed_table = BedTable6Plus(self.extra_field_names, 
                                  extra_column_dtype=self.extra_field_dtype, 
                                  )
        bed_table.load_from_dataframe(self.data_df.copy())
        self.assertArrayEqual(bed_table.to_dataframe().values, self.data_df.values)

    def test_load_from_file(self):
        bed_table = BedTable6Plus(self.extra_field_names,
                                  extra_column_dtype=self.extra_field_dtype,
                                  )
        bed_table.load_from_file(self.data_file)

        self.assertArrayEqual(bed_table.to_dataframe().values, self.data_df.values)

    def __init_test_bed_table(self):
        bed_table = BedTable6Plus(extra_column_names=self.extra_field_names,
                                  extra_column_dtype=self.extra_field_dtype, 
                                  )

        bed_table.load_from_file(self.data_file)

        return bed_table

if __name__ == "__main__":
    unittest.main()
