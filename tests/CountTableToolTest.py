
import argparse
import unittest
import shutil
import sys
import os

import numpy as np
import pandas as pd

sys.path.append("scripts")

from RGTools.BedTable import BedTable6

from scripts.count_table_tool import CountTableTool
from scripts.rsem2count_table import RSEM2CountTable

class CountTableToolTest(unittest.TestCase):
    def setUp(self):
        self.__test_dir = "CountTableToolTest_temp"

        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        # sample data generated from rsem results
        self.__samples = ["sample1", "sample2"]
        self.__rsem_paths = ["sample_data/sample_rsem_results/sample1.rsem.gene.results",
                             "sample_data/sample_rsem_results/sample2.rsem.gene.results",
                             ]
        
        self.__gene_id_list_path = "sample_data/gencode.v37.gene_id.list"

        self.__rsem_count_table = os.path.join(self.__test_dir, "rsem_count_table.csv")
        self.__rsem_count_table_multi_rep = os.path.join(self.__test_dir, "rsem_count_table_multi_rep.csv")
        self.__rsem_region_info = os.path.join(self.__test_dir, "rsem.region_info.bed")

        self._gen_rsem_count_table()

        return super().setUp()

    def tearDown(self):
        if os.path.exists(self.__test_dir):
            shutil.rmtree(self.__test_dir)

        return super().tearDown()

    def _gen_rsem_count_table(self):
        args = argparse.Namespace(sample=self.__samples,
                                  rsem_output=self.__rsem_paths,
                                  gene_id_list=self.__gene_id_list_path,
                                  opath=self.__rsem_count_table, 
                                  count_type="expected_count",
                                  )
        
        RSEM2CountTable.main(args)

        rsem_count_table = CountTableTool.read_input_df(self.__rsem_count_table)

        np.random.seed(76)
        rsem_count_table["sample1_rep1"] = rsem_count_table["sample1"] + np.random.randint(-5, 5, len(rsem_count_table))
        rsem_count_table["sample1_rep2"] = rsem_count_table["sample1"] + np.random.randint(-5, 5, len(rsem_count_table))
        rsem_count_table["sample2_rep1"] = rsem_count_table["sample2"] + np.random.randint(-5, 5, len(rsem_count_table))
        rsem_count_table["sample2_rep2"] = rsem_count_table["sample2"] + np.random.randint(-5, 5, len(rsem_count_table))

        rsem_count_table.drop(columns=["sample1", "sample2"], inplace=True)
        rsem_count_table.clip(lower=0, inplace=True)

        rsem_count_table.to_csv(self.__rsem_count_table_multi_rep)

        # make a temp region info file
        region_info = pd.DataFrame({"chrom": ["chromUn"] * rsem_count_table.shape[0],
                                    "start": range(rsem_count_table.shape[0]), 
                                    "end": range(1, rsem_count_table.shape[0]+1),
                                    "name": ".", 
                                    "score": ".",
                                    "strand": "+",
                                    }, 
                                    index=rsem_count_table.index,
                                    )
        region_info.to_csv(self.__rsem_region_info)

    def get_compute_tissue_tstat_args(self):
        args = argparse.Namespace(subcommand="compute_tissue_tstat",
                                  inpath=self.__rsem_count_table_multi_rep,
                                  tissue_labels="1,1,2,2", 
                                  opath=os.path.join(self.__test_dir, "tissue_tstat.csv"),
                                  )
        
        return args
    
    def test_compute_tissue_tstat_main(self):
        args = self.get_compute_tissue_tstat_args()
        CountTableTool.main(args)

        result_df = CountTableTool.read_input_df(args.opath)
        self.assertAlmostEqual(result_df.loc["ENSG00000185920.16", "1"], -28.831860, places=3)
    
    def get_tstat_table2bed_args(self):
        args = argparse.Namespace(subcommand="tstat_table2bed",
                                  inpath=os.path.join(self.__test_dir, "tissue_tstat.csv"),
                                  percentage=0.2, 
                                  opath=self.__test_dir, 
                                  region_info=self.__rsem_region_info,
                                  )
        
        return args

    def test_tstat_table2bed(self):
        compute_tissue_tstat_args = self.get_compute_tissue_tstat_args()
        CountTableTool.main(compute_tissue_tstat_args)

        tstat_table2bed_args = self.get_tstat_table2bed_args()
        CountTableTool.main(tstat_table2bed_args)

        output_bt = BedTable6()
        output_bt.load_from_file(os.path.join(tstat_table2bed_args.opath, 
                                              "1.tstat_top.20%.bed", 
                                              ))

        self.assertEqual(output_bt.get_start_locs()[1], 19)

