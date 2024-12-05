
import unittest
import argparse
import shutil
import sys
import os

from Bio import SeqIO

sys.path.append("scripts")

from scripts.mutagenesis_TRE import MutaGenesisTRE

class MutaGenesisTRETest(unittest.TestCase):
    def setUp(self):
        self.__test_dir = "MutaGenesisTRETest_temp"
        self.__polymorphism_path = "sample_data/sample.chr2.snp.bed"
        self.__tre_path = "sample_data/sample.snp_overlapping_tre.chr2.bed"
        self.__hg38_path = "large_sample_data/hg38.fa"

        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        return super().setUp()
    
    def tearDown(self):
        shutil.rmtree(self.__test_dir)

        return super().tearDown()

    def __get_default_args(self):
        return argparse.Namespace(inpath_polymorphisms=self.__polymorphism_path, 
                                  inpath_TREs=self.__tre_path, 
                                  opath=os.path.join(self.__test_dir, "out.fa"), 
                                  job_name="test_mutagenesis_tre", 
                                  genome_path=self.__hg38_path, 
                                  TRE_file_type="bed3", 
                                  )
    
    def test_get_sequence(self):
        seq = MutaGenesisTRE.get_sequence(self.__hg38_path, 
                                          "chr2", 
                                          100, 
                                          200, 
                                          )

        self.assertEqual(seq, "N"*100)

        seq = MutaGenesisTRE.get_sequence(self.__hg38_path, 
                                          "chr2", 
                                          100110, 
                                          100120, 
                                          )

        self.assertEqual(seq, "GTGAAAACAC")


    def test_main(self):
        args = self.__get_default_args()

        MutaGenesisTRE.main(args)

        with open(args.opath, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                record_id = record.id
                record_seq = str(record.seq)

                chrom, start, end, variant_info = record_id.split("_")

                if variant_info == "ref":
                    self.assertEqual(record_seq, MutaGenesisTRE.get_sequence(self.__hg38_path, 
                                                                             chrom, 
                                                                             int(start), 
                                                                             int(end), 
                                                                             ))
                
                else:
                    mut_coord = int(variant_info.split(":")[0])
                    mut_bases = variant_info.split(":")[1]
                    ref_base = mut_bases.split("2")[0].upper()
                    alt_base = mut_bases.split("2")[1].upper()

                    self.assertEqual(ref_base, 
                                     MutaGenesisTRE.get_sequence(self.__hg38_path,
                                                                 chrom, 
                                                                 mut_coord, 
                                                                 mut_coord+1, 
                                                                 ).upper())
                    self.assertEqual(alt_base, record_seq[mut_coord-int(start)].upper())