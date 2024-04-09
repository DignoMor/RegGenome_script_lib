
# bed table classes to store bed files as pd.DataFrame

import numpy as np
import pandas as pd


class BedTable3:
    def __init__(self):
        pass
        
    def load_from_file(self, ipath: str) -> None:
        '''
        Load a bed file.
        '''
        pass

    def apply_logical_filter(self, logical_array: np.array) -> BedTable3:
        '''
        Use a logical np.array to filter the table.
        Return a new BedTable3 instance.
        '''
        pass

    def region_subset(self, chrom: str, start: int, end: int) -> BedTable3:
        '''
        Subset the table to the given region.
        Return a new BedTable3 instance.
        '''
        pass

    def to_dataframe(self) -> pd.DataFrame:
        '''
        Return the table as a copy of pd.DataFrame
        '''
        pass

    def write(self, opath: str) -> None:
        '''
        Write the table to a bed file.
        '''
        pass

    def get_chrom_names(self) -> np.array:
        '''
        Return a np.array of chrom names.
        '''
        pass

    def get_start_locs(self) -> np.array:
        '''
        Return a np.array of start locations.
        '''
        pass

    def get_end_locs(self) -> np.array:
        '''
        Return a np.array of end locations.
        '''
        pass
    
    def iter_regions(self) -> tuple:
        '''
        Return a iterator of regions.
        '''
        pass

    def sort(self) -> None:
        '''
        Sort the bed table.
        '''
        pass

    def _is_sort(self) -> bool:
        '''
        Check if the table is sorted.
        '''
        pass
