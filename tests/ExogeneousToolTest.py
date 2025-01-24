
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
        self.__sample_exogeneous_fasta_path = "sample_data/sample_exogeneous_sequence/sample.exogeneous.chr2.fa"
        self.__sample_bed3_path = "sample_data/sample_exogeneous_sequence/sample.exogeneous.chr2.predicted_peaks.bed"

        return super().setUp()
    
    def tearDown(self):
        if os.path.exists(self.__test_dir):
            shutil.rmtree(self.__test_dir)
    
    def get_parse_default_args(self):
        args = argparse.Namespace(fasta=self.__sample_exogeneous_fasta_path,
                                  region_file_path=self.__sample_bed3_path,
                                  region_file_type="bed3",
                                  fasta_out=os.path.join(self.__test_dir, "test_output.fasta"),
                                  region_out=os.path.join(self.__test_dir, "test_output.bed"),
                                  regex="chr2_127105185_127107299.*", 
                                  subcommand="filter",
                                  )
        return args
    
    def test_load_fasta(self):
        seq_ids, seqs = ExogeneousTool.load_fasta(self.__sample_exogeneous_fasta_path)
        
        self.assertEqual(len(seq_ids), 58)
        self.assertEqual(len(seqs), 58)
        self.assertEqual(seq_ids[1], "chr2_127084012_127086126_127084192:C2T")
        self.assertEqual(seqs[1][:5], "ATCAC")

    def test_filter_main(self):
        args = self.get_parse_default_args()
        
        ExogeneousTool.filter_main(args)

        with open(args.fasta_out, "r") as output_f:
            for record in SeqIO.parse(output_f, "fasta"):
                self.assertEqual(record.id[:24], "chr2_127105185_127107299")
                self.assertEqual(record.seq[:5], "CAAAG")
