
import argparse

class IdentifyPausingRegion:
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--tss_bed", 
                            help="Input path for TSS reference file.", 
                            type=str, 
                            required=True, 
                            )
        
        parser.add_argument("--l_pad", 
                            help="Left padding for pausing region. "
                                 "Padded region will be the search space for pausing region.", 
                            type=int, 
                            default=500,
                            )

        parser.add_argument("--r_pad", 
                            help="Right padding for pausing region. "
                                 "Padded region will be the search space for pausing region.", 
                            type=int, 
                            default=499,
                            )
        
        parser.add_argument("--target_size", 
                            help="Target size for pausing region. "
                                 "Pausing region will be the target size region with the highest signal density.", 
                            type=int, 
                            default=250,
                            )
        
        parser.add_argument("--bw_pl", 
                            help="Input path for plus strand bigwig file for PROseq assay.", 
                            type=str,
                            required=True,
                            )

        parser.add_argument("--bw_mn", 
                            help="Input path for minus strand bigwig file for PROseq assay.", 
                            type=str,
                            required=True,
                            )
        
        parser.add_argument("--opath",
                            help="Output path for pausing region bed file.",
                            type=str,
                            required=True,
                            )

    @staticmethod
    def main(args):
        pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify pausing region")

    IdentifyPausingRegion.set_parser(parser)

    args = parser.parse_args()

    IdentifyPausingRegion.main(args)
