
import argparse
import unittest
import shutil
import json
import sys
import os

import numpy as np

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
        self.__sample_pl_track_npy_path = "sample_data/sample_exogeneous_sequence/sample.exogeneous.chr2.pred_profiles_pl.npy"
        self.__sample_mn_track_npy_path = "sample_data/sample_exogeneous_sequence/sample.exogeneous.chr2.pred_profiles_mn.npy"

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

    def get_compare_mutagenesis_default_args(self):
        args = argparse.Namespace(fasta=self.__sample_exogeneous_fasta_path,
                                  region_file_path=self.__sample_bed3_path,
                                  region_file_type="bed3",
                                  regex_ref_seqs="(.*)_ref", 
                                  regex_mut_seqs="(.*)_[0-9]*:[ATCG]2[ATCG]",
                                  pl_track_npy=self.__sample_pl_track_npy_path, 
                                  mn_track_npy=self.__sample_mn_track_npy_path,
                                  opath=self.__test_dir, 
                                  total_count_plot_path=os.path.join(self.__test_dir, "total_count_plot.png"),
                                  )
        
        return args
    
    def test_compare_mutagenesis_main(self):
        args = self.get_compare_mutagenesis_default_args()
        
        ExogeneousTool.compare_mutagenesis_main(args)

        with open(os.path.join(args.opath, "chr2_127104997_127107111.json"), "r") as output_f:
            ref_mut_dict = json.load(output_f)
            self.assertEqual(ref_mut_dict["ref_log_total_count"], 2.607356071472168)
            self.assertEqual(ref_mut_dict["mut_log_total_count"], 
                             [2.607358455657959, 
                              2.6025302410125732, 
                              2.6089253425598145, 
                              2.606109619140625, 
                              ], 
                             )

    def set_up_compute_track_correlation_test(self):
        self.__pseudosample_pl_track_npy_path = os.path.join(self.__test_dir, "pseudosample_pl_track.npy")
        self.__pseudosample_mn_track_npy_path = os.path.join(self.__test_dir, "pseudosample_mn_track.npy")
        pl_track = np.load(self.__sample_pl_track_npy_path)
        mn_track = np.load(self.__sample_mn_track_npy_path)

        np.random.seed(76)
        pl_track_rand = pl_track + np.random.normal(0, 0.1, pl_track.shape)
        mn_track_rand = mn_track + np.random.normal(0, 0.1, mn_track.shape)

        pl_track_rand[pl_track_rand < 0] = 0
        mn_track_rand[mn_track_rand > 0] = 0

        np.save(self.__pseudosample_pl_track_npy_path, pl_track_rand)
        np.save(self.__pseudosample_mn_track_npy_path, mn_track_rand)

    def get_compute_track_correlation_default_args(self):
        args = argparse.Namespace(
            subcommand="compute_track_correlation",
            pl_track1_npy=self.__sample_pl_track_npy_path,
            mn_track1_npy=self.__sample_mn_track_npy_path,
            pl_track2_npy=self.__pseudosample_pl_track_npy_path,
            mn_track2_npy=self.__pseudosample_mn_track_npy_path,
            jensenshannon=os.path.join(self.__test_dir, "js_dist.npy"),
        )
        
        return args

    def test_compute_track_correlation(self):
        self.set_up_compute_track_correlation_test()
        args = self.get_compute_track_correlation_default_args()
        ExogeneousTool.main(args)

        js_dist = np.load(args.jensenshannon)
        self.assertEqual(js_dist.shape, (116, ))
        self.assertAlmostEqual(js_dist[0], 0.281159880, places=4)
