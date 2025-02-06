
import argparse
import unittest
import shutil
import sys
import os

import numpy as np

sys.path.append("scripts")

from scripts.count_table_tool import CountTableTool
from scripts.rsem2count_table import RSEM2CountTable

class CountTableToolTest(unittest.TestCase):
    def setUp(self):
        self.__test_dir = "CountTableToolTest_temp"

        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        self.__samples = ["sample1", "sample2"]
        self.__rsem_paths = ["sample_data/sample_rsem_results/sample1.rsem.gene.results",
                             "sample_data/sample_rsem_results/sample2.rsem.gene.results",
                             ]
        
        self.__gene_id_list_path = "sample_data/gencode.v37.gene_id.list"

        self.__sample_count_table = os.path.join(self.__test_dir, "sample_count_table.tsv")
        self.__sample_count_table_multi_rep = os.path.join(self.__test_dir, "sample_count_table_multi_rep.tsv")

        self._gen_sample_count_table()

        return super().setUp()

    def tearDown(self):
        if os.path.exists(self.__test_dir):
            shutil.rmtree(self.__test_dir)

        return super().tearDown()

    def _gen_sample_count_table(self):
        args = argparse.Namespace(sample=self.__samples,
                                  rsem_output=self.__rsem_paths,
                                  gene_id_list=self.__gene_id_list_path,
                                  opath=self.__sample_count_table, 
                                  count_type="expected_count",
                                  )
        
        RSEM2CountTable.main(args)

        sample_count_table = CountTableTool.read_input_df(self.__sample_count_table)

        np.random.seed(76)
        sample_count_table["sample1_rep1"] = sample_count_table["sample1"] + np.random.randint(-5, 5, len(sample_count_table))
        sample_count_table["sample1_rep2"] = sample_count_table["sample1"] + np.random.randint(-5, 5, len(sample_count_table))
        sample_count_table["sample2_rep1"] = sample_count_table["sample2"] + np.random.randint(-5, 5, len(sample_count_table))
        sample_count_table["sample2_rep2"] = sample_count_table["sample2"] + np.random.randint(-5, 5, len(sample_count_table))

        sample_count_table.drop(columns=["sample1", "sample2"], inplace=True)
        sample_count_table.clip(lower=0, inplace=True)

        sample_count_table.to_csv(self.__sample_count_table_multi_rep)

    def get_compute_tissue_tstat_args(self):
        args = argparse.Namespace(subcommand="compute_tissue_tstat",
                                  inpath=self.__sample_count_table_multi_rep,
                                  tissue_labels="1,1,2,2", 
                                  opath=os.path.join(self.__test_dir, "tissue_tstat.tsv"),
                                  )
        
        return args
    
    def test_compute_tissue_tstat_main(self):
        args = self.get_compute_tissue_tstat_args()
        CountTableTool.main(args)

        result_df = CountTableTool.read_input_df(args.opath)
        self.assertAlmostEqual(result_df.loc["ENSG00000185920.16", "1"], -28.831860, places=3)
