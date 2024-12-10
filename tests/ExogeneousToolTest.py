
import argparse
import unittest
import shutil
import sys
import os

from Bio import SeqIO

sys.path.append("scripts")
from scripts.exogeneous_tool import ExogeneousTool

class ExogeneousToolTest(unittest.TestCase):
    def setUp(self):
        self.__test_dir = "ExogeneousToolTest_temp"

        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        self.__sample_fasta_path = "large_sample_data/hg38.fa"

        return super().setUp()
    
    def tearDown(self):
        if os.path.exists(self.__test_dir):
            shutil.rmtree(self.__test_dir)
    
    def get_parse_default_args(self):
        args = argparse.Namespace(inpath=self.__sample_fasta_path,
                                  outpath=os.path.join(self.__test_dir, "test_output.fasta"),
                                  regex="chr11_KI270721v1_*", 
                                  subcommand="filter",
                                  )
        return args

    def test_filter_main(self):
        args = self.get_parse_default_args()
        
        ExogeneousTool.filter_main(args)

        with open(args.outpath, "r") as output_f:
            for record in SeqIO.parse(output_f, "fasta"):
                self.assertEqual(record.id, "chr11_KI270721v1_random")
                self.assertEqual(record.seq[:5], "ctgcg")
