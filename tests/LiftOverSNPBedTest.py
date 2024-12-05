
import argparse
import unittest
import shutil
import sys
import os

sys.path.append("scripts")

from scripts.lift_over_snp_bed import LiftOverSNPBed

class LiftOverSNPBedTest(unittest.TestCase):
    def setUp(self):
        self.__test_dir = "LiftOverSNPBedTest_temp_data"
        self.__sample_snp_bed = os.path.join("sample_data/sample.chr2.snp.bed")
        self.__short_sample_snp_bed = os.path.join(self.__test_dir, 
                                                   "sample.chr2.10.snp.bed", 
                                                   )

        self.file_head(self.__sample_snp_bed, 
                       self.__short_sample_snp_bed,
                       n=10, 
                       )
        if not os.path.exists(self.__test_dir):
            os.makedirs(self.__test_dir)

        return super().setUp()
    

    def tearDown(self):
        if os.path.exists(self.__test_dir):
            shutil.rmtree(self.__test_dir)
        return super().tearDown()
    
    def file_head(self, inpath, opath, n=10):
        '''
        Write the first n lines of inpath to opath.
        '''
        with open(inpath, "r") as infile:
            with open(opath, "w") as outfile:
                for i in range(n):
                    line = infile.readline()
                    outfile.write(line)
    
    def get_simple_args(self):
        args = argparse.Namespace(log_file=os.path.join(self.__test_dir, "test.log"),
                                  bed_in=self.__short_sample_snp_bed, 
                                  from_genome="GRCh38",
                                  to_genome="hg19",
                                  bed_out=os.path.join(self.__test_dir, "test.out.bed"),
                                  )
        return args
    
    def test_main(self):
        args = self.get_simple_args()

        LiftOverSNPBed.main(args)

        output_bt = LiftOverSNPBed.SNPBedTable()
        output_bt.load_from_file(args.bed_out)

        input_bt = LiftOverSNPBed.SNPBedTable()
        input_bt.load_from_file(args.bed_in)

        self.assertTrue((input_bt.get_chrom_names() == output_bt.get_chrom_names()).all())
        self.assertTrue((input_bt.get_region_names() == output_bt.get_region_names()).all())
        self.assertEqual(output_bt.get_start_locs()[0], 106366055)


    



