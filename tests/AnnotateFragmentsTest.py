
import unittest
import requests
import argparse
import shutil
import gzip
import os

import pandas as pd
import numpy as np

from scripts.annotate_fragments import AnnotateFragments

class AnnotateFragmentsTest(unittest.TestCase):
    def setUp(self):
        self._test_dir = "AnnotateFragmentsTest_temp"
        self.fragment_file = os.path.join(self._test_dir, "fragments.tsv")
        self.annotation_file = os.path.join(self._test_dir, "annotation.tsv")

        if not os.path.exists(self._test_dir):
            os.makedirs(self._test_dir)
        
        self.download_test_data()
    
    def download_test_data(self):
        # fragments
        url = "https://personal.broadinstitute.org/bjames/AD_snATAC/fragments/D19-13166_fragments.bed.gz"
        response = requests.get(url)
        with open(self.fragment_file + ".gz", "wb") as f:
            f.write(response.content)
        with gzip.open(self.fragment_file + ".gz", "rb") as f_in:
            with open(self.fragment_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        # annotation
        url = "https://personal.broadinstitute.org/bjames/AD_snATAC/integration/integration.tsv.gz"
        response = requests.get(url)
        with open(self.annotation_file + ".gz", "wb") as f:
            f.write(response.content)
        with gzip.open(self.annotation_file + ".gz", "rb") as f_in:
            with open(self.annotation_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        anno_df = pd.read_csv(self.annotation_file + ".gz", 
                              sep="\t", 
                              compression="gzip", 
                              index_col=0, 
                              )

        cell_type_anno_df = anno_df[["MajorCellType"]].copy()
        cell_type_anno_df.replace("Opc", "OPC", inplace=True)
        cell_type_anno_df = cell_type_anno_df.loc[np.array([s.startswith("D19-13166") for s in cell_type_anno_df.index])] 
        cell_type_anno_df.index = cell_type_anno_df.index.str.replace("D19-13166#", "")

        cell_type_anno_df.to_csv(self.annotation_file, 
                                 sep="\t", 
                                 header=False, 
                                 )

        os.remove(self.annotation_file + ".gz")
        os.remove(self.fragment_file + ".gz")

    def tearDown(self):
        if os.path.exists(self._test_dir):
            shutil.rmtree(self._test_dir)
        super().tearDown()

    def test_annotate_fragments(self):
        args = argparse.Namespace(
            fragments=self.fragment_file,
            annotation=self.annotation_file,
            oheader=os.path.join(self._test_dir, "test_output"),
        )

        AnnotateFragments.main(args)

        exc_output = pd.read_csv(args.oheader + ".Exc.tsv", 
                                 sep="\t", 
                                 header=None, 
                                 names=["chrom", "start", "end", "cell_id", "count"],
                                 )
        
        self.assertEqual(exc_output.shape[0], 12808521)

        opc_output = pd.read_csv(args.oheader + ".OPC.tsv", 
                                 sep="\t", 
                                 header=None, 
                                 names=["chrom", "start", "end", "cell_id", "count"],
                                 )
        self.assertEqual(opc_output.shape[0], 1245186)
