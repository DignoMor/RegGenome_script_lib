
import unittest
import argparse
import shutil
import sys
import os

import pandas as pd
import numpy as np

from scripts.RGTools.BedTable import BedTable3, BedTable6, BedTable6Plus

sys.path.append('scripts')
from scripts.genomicelement_tool import GenomicElementTool

class GenomicElementToolTest(unittest.TestCase):
    def setUp(self):
        self.__temp_dir = "GenomicElementToolTest_temp_data"

        if not os.path.exists(self.__temp_dir):
            os.makedirs(self.__temp_dir)

        self.setup_count_bw_test()
        
        return super().setUp()

    def setup_count_bw_test(self):
        self.__bw_pls = ["sample_data/ENCFF993VCR.pl.bigWig"]
        self.__bw_mns = ["sample_data/ENCFF182TPF.mn.bigWig"]

        self.__bed3_path = os.path.join(self.__temp_dir, "test.bed3")
        self.__bed6_path = os.path.join(self.__temp_dir, "test.bed6")
        self.__bed6gene_path = os.path.join(self.__temp_dir, "test.bed6gene")

        region_df = pd.DataFrame({"chr": ["chr6", "chr14", "chr17"],
                                  "start_loc": [170553801, 75278325, 45894026],
                                  "end_loc": [170554802, 75279326, 45895027],
                                  "name": ["TBP", "FOS", "MAPT"],
                                  "score": [1, 2, 3],
                                  "strand": ["+", "-", "+"],
                                  "gene_name": ["gene1", "gene2", "gene3"],
                                  }, 
                                  )

        # Write regions in different bed format
        bedtable3 = BedTable3()
        bedtable3.load_from_dataframe(region_df[["chr", "start_loc", "end_loc"]], 
                                      column_map={"chrom": "chr", 
                                                  "start": "start_loc", 
                                                  "end": "end_loc", 
                                                  },
                                      )
        bedtable3.write(self.__bed3_path)

        bedtable6 = BedTable6()
        bedtable6.load_from_dataframe(region_df[["chr", "start_loc", "end_loc", "name", "score", "strand"]], 
                                      column_map={"chrom": "chr", 
                                                  "start": "start_loc", 
                                                  "end": "end_loc", 
                                                  "name": "name", 
                                                  "score": "score", 
                                                  "strand": "strand", 
                                                  },
                                      )
        bedtable6.write(self.__bed6_path)

        bedtable6gene = BedTable6Plus(extra_column_names=["gene_name"],
                                      extra_column_dtype=[str],
                                      )
        bedtable6gene.load_from_dataframe(region_df, 
                                          column_map={"chrom": "chr", 
                                                      "start": "start_loc", 
                                                      "end": "end_loc", 
                                                      "name": "name", 
                                                      "score": "score", 
                                                      "strand": "strand", 
                                                      "gene_name": "gene_name", 
                                                      },
                                          )
        bedtable6gene.write(self.__bed6gene_path)

    def tearDown(self):
        if os.path.exists(self.__temp_dir):
            shutil.rmtree(self.__temp_dir)

        return super().tearDown()
    
    def get_count_bw_simple_args(self):
        args = argparse.Namespace()
        args.subcommand = "count_bw"
        args.bw_pl = self.__bw_pls[0]
        args.bw_mn = self.__bw_mns[0]
        args.single_bw = False
        args.region_file_path = self.__bed6_path
        args.region_file_type = "bed6"
        args.genome_path = None
        args.quantification_type = "raw_count"
        args.opath = os.path.join(self.__temp_dir, "output.npy")

        return args

    def test_count_bw(self):
        args = self.get_count_bw_simple_args()

        GenomicElementTool.count_bw_main(args)

        output = np.load(args.opath)

        self.assertEqual(output.shape, (3,))
        self.assertEqual(output[0], 17)
        self.assertEqual(output[2], 348)


