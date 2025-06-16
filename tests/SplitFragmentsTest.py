
import os
import gzip
import shutil
import argparse
import unittest
import requests

import pandas as pd

from scripts.split_fragments import SplitFragments


class SplitFragmentsTest(unittest.TestCase):
    def setUp(self):
        self.__test_dir = "SplitFragmentsTest"
        os.makedirs(self.__test_dir, exist_ok=True)

        self.__fragments_file = os.path.join(self.__test_dir, "fragments.tsv.gz")
        self.__label_file = os.path.join(self.__test_dir, "labels.tsv")

        self.__download_example_data()
    
    def __download_example_data(self):
        url = "https://personal.broadinstitute.org/bjames/AD_snATAC/fragments/D19-8612_fragments.bed.gz"
        response = requests.get(url)
        with open(self.__fragments_file, "wb") as f:
            f.write(response.content)
        
        url = "https://personal.broadinstitute.org/bjames/AD_snATAC/integration/integration.tsv.gz"
        response = requests.get(url)
        with open(self.__label_file + ".temp.tsv.gz", "wb") as f:
            f.write(response.content)

        label_df = pd.read_csv(self.__label_file + 
                               ".temp.tsv.gz", 
                               sep="\t", 
                               compression="gzip",
                               )
        label_df = label_df.iloc[:, [0, 6]]
        label_df.iloc[:, 0] = label_df.iloc[:, 0].agg(lambda x: x.split("#")[-1])
        label_df.iloc[:, 1] = label_df.iloc[:, 1].replace("Opc", "OPC")
        label_df.to_csv(self.__label_file, 
                        sep="\t", 
                        index=False,
                        header=False,
                        )

        os.remove(self.__label_file + ".temp.tsv.gz")


    def tearDown(self):
        shutil.rmtree(self.__test_dir)
        super().tearDown()

    def test_split_fragments(self):
        args = argparse.Namespace(
            fragments=self.__fragments_file,
            label_file=self.__label_file,
            oheader=os.path.join(self.__test_dir, "split_fragments"),
        )

        SplitFragments.main(args)

        exc_df = pd.read_csv(os.path.join(self.__test_dir, "split_fragments.Exc.fragments.tsv"), 
                             sep="\t", 
                             header=None, 
                             )
        self.assertEqual(exc_df.shape[0], 1669846)
        self.assertEqual(exc_df.shape[1], 5)

if __name__ == "__main__":
    unittest.main()
