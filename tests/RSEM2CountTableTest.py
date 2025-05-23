
import unittest
import argparse
import shutil
import sys
import os

import pandas as pd

sys.path.append("scripts")
from scripts.rsem2count_table import RSEM2CountTable

class RSEM2CountTableTest(unittest.TestCase):
    def setUp(self):
        self.__wdir = "RSEM2CountTableTest_temp_data"

        self.__samples = ["sample1", "sample2"]
        self.__rsem_paths = ["sample_data/sample_rsem_results/sample1.rsem.gene.results",
                             "sample_data/sample_rsem_results/sample2.rsem.gene.results",
                             ]
        
        self.__gene_id_list_path = "sample_data/gencode.v37.gene_id.list"

        if not os.path.exists(self.__wdir):
            os.makedirs(self.__wdir)

        return super().setUp()
    
    def tearDown(self):
        shutil.rmtree(self.__wdir)

        return super().tearDown()
    
    def get_simple_args(self):
        args = argparse.Namespace(sample=self.__samples,
                                  rsem_output=self.__rsem_paths,
                                  gene_id_list=self.__gene_id_list_path,
                                  opath=os.path.join(self.__wdir, "count_table.tsv"),
                                  count_type="expected_count",
                                  )
        
        return args
    
    def test_main(self):
        args = self.get_simple_args()
        
        RSEM2CountTable.main(args)

        output = pd.read_csv(args.opath, 
                             index_col=0, 
                             )
        
        self.assertEqual(output.loc["ENSG00000243485.5", "sample1"], 2.69)
    
    def test_tpm_output(self):
        args = self.get_simple_args()
        args.count_type = "TPM"
        
        RSEM2CountTable.main(args)

        output = pd.read_csv(args.opath, 
                             index_col=0, 
                             )
        
        self.assertEqual(output.loc["ENSG00000243485.5", "sample1"], 0.4)

