
# bed table classes to store bed files as pd.DataFrame

import numpy as np
import pandas as pd

from .exceptions import BedTableLoadException

class BedTableIterator:
    def __init__(self, bed_table: 'BedTable3') -> None:
        self.__bed_table = bed_table
        self.__index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.__index >= len(self.__bed_table):
            raise StopIteration
        else:
            region = self.__bed_table.get_region_by_index(self.__index)
            self.__index += 1
            return region

class BedTable3:
    def __init__(self):
        self.__column_names = ["chrom", "start", "end"]
        self.__data_df = pd.DataFrame(columns=self.__column_names)
        
    def load_from_file(self, ipath: str) -> None:
        '''
        Load a bed file.
        '''
        try:
            self.__data_df = pd.read_csv(ipath, 
                                         sep="\t", 
                                         names=self.__column_names,
                                         )
        except ValueError as e:
            raise BedTableLoadException(f"Error loading bed file: number of columns does not match.")
        
        if not self._is_sorted():
            self._sort()

    def load_from_dataframe(self, df: pd.DataFrame) -> None:
        '''
        Load a pd.DataFrame.
        '''
        try:
            self.__data_df = pd.DataFrame(df.values, columns=self.__column_names)
        except ValueError as e:
            raise BedTableLoadException(f"Error loading pd.DataFrame: number of columns does not match.")

        if not self._is_sorted():
            self._sort()

    def apply_logical_filter(self, logical_array: np.array) -> 'BedTable3':
        '''
        Use a logical np.array to filter the table.
        Return a new BedTable3 instance.
        '''
        if not isinstance(logical_array, np.ndarray):
            raise ValueError("logical_array must be a np.array")
        if not logical_array.dtype == bool:
            raise ValueError("logical_array must be a boolean np.array")

        new_bed_table = BedTable3()
        new_bed_table.load_from_dataframe(self.__data_df.loc[logical_array])

        return new_bed_table
        

    def region_subset(self, chrom: str, start: int, end: int) -> 'BedTable3':
        '''
        Subset the table to the given region.
        Only return regions that are fully contained in the given region.
        Return a new BedTable3 instance.
        '''
        subset_data_df = self.__data_df.loc[self.__data_df["chrom"] == chrom]
        subset_data_df = subset_data_df.loc[subset_data_df["start"] >= start]
        subset_data_df = subset_data_df.loc[subset_data_df["end"] <= end]

        new_bed_table = BedTable3()
        new_bed_table.load_from_dataframe(subset_data_df)

        return new_bed_table

    def to_dataframe(self) -> pd.DataFrame:
        '''
        Return the table as a copy of pd.DataFrame
        '''
        return self.__data_df.copy()

    def write(self, opath: str) -> None:
        '''
        Write the table to a bed file.
        '''
        self.__data_df.to_csv(opath, 
                              sep="\t", 
                              header=False, 
                              index=False, 
                              )

    def get_chrom_names(self) -> np.array:
        '''
        Return a np.array of chrom names.
        '''
        return self.__data_df["chrom"].values

    def get_start_locs(self) -> np.array:
        '''
        Return a np.array of start locations.
        '''
        return self.__data_df["start"].values

    def get_end_locs(self) -> np.array:
        '''
        Return a np.array of end locations.
        '''
        return self.__data_df["end"].values
    
    def get_region_by_index(self, index: int) -> pd.Series:
        '''
        Return a region by index.
        '''
        return self.__data_df.iloc[index]
    
    def iter_regions(self) -> tuple:
        '''
        Return a iterator of regions.
        '''
        return BedTableIterator(self)

    def _sort(self) -> None:
        '''
        Sort the bed table.
        '''
        self.__data_df.sort_values(by=["chrom", "start", "end"], inplace=True)
        self.__data_df.reset_index(drop=True, inplace=True)

    def _is_sorted(self) -> bool:
        '''
        Check if the table is sorted.
        '''
        for i in range(1, len(self.__data_df)):
            if self.__data_df.iloc[i]["chrom"] < self.__data_df.iloc[i-1]["chrom"]:
                return False
            if self.__data_df.iloc[i]["chrom"] == self.__data_df.iloc[i-1]["chrom"]:
                if self.__data_df.iloc[i]["start"] < self.__data_df.iloc[i-1]["start"]:
                    return False
                if self.__data_df.iloc[i]["start"] == self.__data_df.iloc[i-1]["start"]:
                    if self.__data_df.iloc[i]["end"] < self.__data_df.iloc[i-1]["end"]:
                        return False

        return True

    def __len__(self) -> int:
        return self.__data_df.shape[0]