
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
        super().__init__()
        self._data_df = pd.DataFrame(columns=self.column_names)

    # public methods
    @property
    def column_names(self):
        return ["chrom", "start", "end"]
    
    @property
    def column_types(self):
        '''
        A dictionary of column types.
        '''
        return {"chrom": str, "start": int, "end": int}

    def load_from_file(self, ipath: str) -> None:
        '''
        Load a bed file.
        '''
        try:
            self._data_df = pd.read_csv(ipath, 
                                         sep="\t", 
                                         names=self.column_names,
                                         )
        except ValueError as e:
            raise BedTableLoadException(f"Error loading bed file: number of columns does not match.")
        
        self.__force_dtype()

        if not self._is_sorted():
            self._sort()

    def load_from_dataframe(self, df: pd.DataFrame, 
                            column_map=None) -> None:
        '''
        Load a pd.DataFrame.

        Keyword arguments:
        df -- a pd.DataFrame that contains the data for a bed table
        column_map -- a dictionary that maps the column names in the 
                      bed table to the column names in the pd.DataFrame         
        '''
        if not column_map:
            column_map = {col: col for col in self.column_names}

        try:
            self._data_df = pd.DataFrame(df[[column_map[col] for col in self.column_names]].values, 
                                         columns=self.column_names, 
                                         )
        except ValueError as e:
            raise BedTableLoadException(f"Error loading pd.DataFrame: number of columns does not match.")
        
        self.__force_dtype()

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

        new_bed_table = self.__class__()
        new_bed_table.load_from_dataframe(self._data_df.loc[logical_array])

        return new_bed_table
        
    def region_subset(self, chrom: str, start: int, end: int) -> 'BedTable3':
        '''
        Subset the table to the given region.
        Only return regions that are fully contained in the given region.
        Return a new BedTable3 instance.
        '''
        #TODO: Add only subsetting by chrome option
        subset_data_df = self._data_df.loc[self._data_df["chrom"] == chrom]
        subset_data_df = subset_data_df.loc[subset_data_df["start"] >= start]
        subset_data_df = subset_data_df.loc[subset_data_df["end"] <= end]

        new_bed_table = self.__class__()
        new_bed_table.load_from_dataframe(subset_data_df)

        return new_bed_table

    def to_dataframe(self) -> pd.DataFrame:
        '''
        Return the table as a copy of pd.DataFrame
        '''
        return self._data_df.copy()

    def write(self, opath: str) -> None:
        '''
        Write the table to a bed file.
        '''
        self._data_df.to_csv(opath, 
                              sep="\t", 
                              header=False, 
                              index=False, 
                              )

    def get_chrom_names(self) -> np.array:
        '''
        Return a np.array of chrom names.
        '''
        return self._data_df["chrom"].values

    def get_start_locs(self) -> np.array:
        '''
        Return a np.array of start locations.
        '''
        return self._data_df["start"].values

    def get_end_locs(self) -> np.array:
        '''
        Return a np.array of end locations.
        '''
        return self._data_df["end"].values
    
    def get_region_by_index(self, index: int) -> pd.Series:
        '''
        Return a region by index.
        '''
        return self._data_df.iloc[index]
    
    def iter_regions(self) -> tuple:
        '''
        Return a iterator of regions.
        '''
        return BedTableIterator(self)

    def _sort(self) -> None:
        '''
        Sort the bed table.
        '''
        self._data_df.sort_values(by=["chrom", "start", "end"], inplace=True)
        self._data_df.reset_index(drop=True, inplace=True)

    def _is_sorted(self) -> bool:
        '''
        Check if the table is sorted.
        '''
        for i in range(1, len(self._data_df)):
            if self._data_df.iloc[i]["chrom"] < self._data_df.iloc[i-1]["chrom"]:
                return False
            if self._data_df.iloc[i]["chrom"] == self._data_df.iloc[i-1]["chrom"]:
                if self._data_df.iloc[i]["start"] < self._data_df.iloc[i-1]["start"]:
                    return False
                if self._data_df.iloc[i]["start"] == self._data_df.iloc[i-1]["start"]:
                    if self._data_df.iloc[i]["end"] < self._data_df.iloc[i-1]["end"]:
                        return False

        return True

    def __force_dtype(self):
        '''
        Force the column types.
        '''
        for field, field_dtype in self.column_types.items():
            self._data_df[field] = self._data_df[field].astype(field_dtype)

    def __len__(self) -> int:
        return self._data_df.shape[0]

class BedTable6(BedTable3):
    def __init__(self):
        super().__init__()

    @property
    def column_names(self):
        return ["chrom", "start", "end", "name", "score", "strand"]

    @property
    def column_types(self):
        column_type = super().column_types
        column_type["name"] = str
        column_type["score"] = float
        column_type["strand"] = str
        return column_type
    
    def get_region_names(self) -> np.array:
        '''
        Return a np.array of region names.
        '''
        return self._data_df["name"].values
    
    def get_region_scores(self) -> np.array:
        '''
        Return a np.array of region scores.
        '''
        return self._data_df["score"].values
    
    def get_region_strands(self) -> np.array:
        '''
        Return a np.array of region strands.
        '''
        return self._data_df["strand"].values

class BedTable6Plus(BedTable6):
    def __init__(self, extra_column_names: list, extra_column_dtype=None):
        self._extra_column_names = extra_column_names

        if not extra_column_dtype:
            self._extra_column_dtype = [str] * len(extra_column_names)
        else:
            self._extra_column_dtype = extra_column_dtype

        super().__init__()

    @property
    def column_names(self):
        return ["chrom", "start", "end", "name", "score", "strand"] + self._extra_column_names
    
    @property
    def column_types(self):
        column_type = super().column_types
        for extra_col, extra_col_dtype in zip(self._extra_column_names, self._extra_column_dtype):
            column_type[extra_col] = extra_col_dtype
        
        return column_type

    def get_region_extra_column(self, column_name) -> np.array:
        '''
        Return a np.array of extra column data for all the regions. Given the column name.
        '''
        return self._data_df[column_name].values
