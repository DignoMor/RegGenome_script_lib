#!/usr/bin/env python3

import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from RGTools.BedTable import BedTable3, BedTable6, BedTable6Plus

class MutaGenesisTRE:
    @staticmethod
    def set_parser(parser):
        parser.add_argument("--inpath_polymorphisms", 
                            help="Input bed6+ file for polymorphisms.", 
                            required=True, 
                            )
        
        parser.add_argument("--inpath_TREs",
                            help="Input bed file for transcription regulatory elements.",
                            required=True,
                            )
        
        parser.add_argument("--opath",
                            help="Output path for mutated TREs.",
                            required=True,
                            )
        
        parser.add_argument("--job_name", 
                            help="Name of the job.",
                            required=True,
                            )
        
        parser.add_argument("--genome_path", 
                            help="Path to the genome fasta file.",
                            required=True,
                            )
        
        parser.add_argument("--TRE_file_type", 
                            help="Type of the TRE file. [bed3] (Choices: {})".format(
                                ", ".join(MutaGenesisTRE.get_supported_TRE_filetype2class_dict().keys())
                            ),
                            default="bed3",
                            )

    @staticmethod
    def get_supported_TRE_filetype2class_dict():
        return {"bed3": BedTable3, 
                "bed6": BedTable6,
                }
    
    @staticmethod
    def input_check(args):
        if not args.TRE_file_type in MutaGenesisTRE.get_supported_TRE_filetype2class_dict():
            raise ValueError("Invalid TRE file type: {}".format(args.TRE_file_type))
    
    @staticmethod
    def load_TRE(tre_path, file_type):
        input_bt = MutaGenesisTRE.get_supported_TRE_filetype2class_dict()[file_type]()

        input_bt.load_from_file(tre_path)
        return input_bt
    
    @staticmethod
    def BedTablePolymorphism():
        '''
        Constructor for BedTablePolymorphism class.
        '''
        bt = BedTable6Plus(extra_column_names=["bases"], 
                           extra_column_dtype=[str], 
                           )
        
        return bt
    
    @staticmethod
    def get_sequence(fasta_path, chrom, start, end):

        with open(fasta_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if record.id == chrom:
                    # Extract the sequence (convert start/end to 0-based indexing)
                    seq = str(record.seq[start:end])
                    break

        return seq
    
    @staticmethod
    def mutagenesis(ref_sequence: str, tre_region: "BedRegion", 
                    index2mut: int, mutated_base: str):
        '''
        Mutagenesis on a ref_sequence and return the mutated sequence.
        
        Keyword arguments:
        - ref_sequence: Reference sequence.
        - tre_region: TRE BedRegion.
        - index2mut: Index of the string to mutate.
        - mutated_base: Base to mutate to

        Return the mutated sequence.
        '''
        mutated_seq = ref_sequence[:index2mut] + mutated_base + ref_sequence[index2mut+1:]
        return mutated_seq
    
    @staticmethod
    def main(args):
        MutaGenesisTRE.input_check(args)

        polymorphisms_bt = MutaGenesisTRE.BedTablePolymorphism()
        polymorphisms_bt.load_from_file(args.inpath_polymorphisms)

        TREs_bt = MutaGenesisTRE.load_TRE(args.inpath_TREs, 
                                          args.TRE_file_type, 
                                          )
        
        if os.path.exists(args.opath):
            os.remove(args.opath)
        
        for tre_region in TREs_bt.iter_regions():

            overlapping_snp_bt = polymorphisms_bt.region_subset(tre_region["chrom"],
                                                                tre_region["start"],
                                                                tre_region["end"],
                                                                )
            if len(overlapping_snp_bt) == 0:
                continue

            ref_sequence = MutaGenesisTRE.get_sequence(args.genome_path, 
                                                       tre_region["chrom"], 
                                                       tre_region["start"], 
                                                       tre_region["end"], 
                                                       )
            
            ref_seq_id = "_".join([
                tre_region["chrom"],
                str(tre_region["start"]),
                str(tre_region["end"]),
                "ref", 
            ])

            record = SeqRecord(Seq(ref_sequence), id=ref_seq_id, description="")
        
            with open(args.opath, "a") as fasta_out:
                SeqIO.write(record, fasta_out, "fasta")

            for snp in overlapping_snp_bt.iter_regions():
                
                index2mut = snp["start"] - tre_region["start"]
                ref_base = ref_sequence[index2mut].upper()
                for base in snp["bases"].split("/"):
                    mutated_base = base.upper()
                    if mutated_base == ref_base:
                        continue
                    if len(mutated_base) > 1:
                        continue

                    mutated_seqs = MutaGenesisTRE.mutagenesis(ref_sequence,
                                                              tre_region, 
                                                              index2mut, 
                                                              mutated_base, 
                                                              )
                    
                    mutated_seqs_id = "_".join([
                        tre_region["chrom"],
                        str(tre_region["start"]).upper(),
                        str(tre_region["end"]).upper(),
                        "{:d}:{}2{}".format(snp["start"], ref_base, mutated_base), 
                    ])

                    record = SeqRecord(Seq(mutated_seqs), id=mutated_seqs_id, description="")
                    with open(args.opath, "a") as fasta_out:
                        SeqIO.write(record, fasta_out, "fasta")
                

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mutagenesis of Transcription Regulatory Elements")
    MutaGenesisTRE.set_parser(parser)
    args = parser.parse_args()

    MutaGenesisTRE.main(args)

