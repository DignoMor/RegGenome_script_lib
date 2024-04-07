#!/usr/bin/env python

# Edit a fasta file and produce a new fasta

import sys
import argparse

import pandas as pd
import numpy as np

from Bio import SeqIO, Seq, SeqRecord

def set_parser(parser):
    parser.add_argument("--ref_fasta", 
                        help="Path to the reference fasta file.", 
                        required=True, 
                        type=str, 
                        )
    
    parser.add_argument("--opath", 
                        help="Path to the output edited fasta file.", 
                        required=True, 
                        type=str, 
                        )
    
    parser.add_argument("--base_info_path",
                        help="Path to the base info file. "
                             "File should be in standard bed6 format with "
                             "an extra 7th column for bases. Bases should be " 
                             "comma seperated. Currently only suport 2 polymorphisms.", 
                        required=True, 
                        type=str, 
                        )
    
    parser.add_argument("--filtered_base_info_opath", 
                        help="Output path for filtered base info file.", 
                        default=None, 
                        type=str, 
                        )

def mutate_seq_seg(seq, start_loc, end_loc, edit_loc, bases):
    '''
    Mutate the base defined by start_loc in the given seq and return a mutated new seq instance.
    The returned new seq would be the segment defined by start_loc and end_loc in 
    python slicing convension.
    Exception will be raised if the reference base is not one of the given bases.

    Keyword Arguments: 
    seq - reference sequence
    start_loc - start location, 0-based.
    end_loc - end location, 0-based.
    start_loc - location for the mutation, 0-based index.
    bases - comma separated list of 2 bases. The reference bases must be one of the 2.
    '''
    ref_base = seq[edit_loc].upper()
    bases = bases.upper().split(",")

    if not len(bases) == 2:
        raise Exception("Only support 2 bases!")
    
    if not ref_base in bases:
        raise Exception("The ref base must be one of the 2 bases!")
    
    bases.remove(ref_base)
    alt_base = bases[0]

    new_seq = seq[start_loc:edit_loc] + alt_base + seq[edit_loc+1:end_loc]

    return new_seq

def mutate_contig(ref_seq, base_info_df):
    '''
    Mutate a whole contig according to the base info df
    '''
    new_seq = Seq.Seq("")
    for _, base_info_row in base_info_df.loc[base_info_df["chrom"] == contig_name].iterrows():
        new_seq += mutate_seq_seg(ref_seq, 
                                  start_loc=len(new_seq), 
                                  end_loc=base_info_row["start_loc"], 
                                  edit_loc=base_info_row["start_loc"] - 1, 
                                  bases=base_info_row["bases"], 
                                  )

    new_seq += ref_seq[len(new_seq):]
    return new_seq

def filter_base_info(base_info_df):
    '''
    Filter the base info df and return a copy df
    '''
    # length of ref and alt base

    # filter for number of polymorphisms
    base_lists = base_info_df["bases"].agg(lambda s: s.split(","))
    num_poly = base_lists.agg(len)
    num_poly_filter = num_poly == 2

    sys.stderr.write("{:d} polymorphisms filtered due to more than 2 alleles.\n".format((~num_poly_filter).sum()))

    # filter insersion/deletions
    base_info_df = base_info_df.loc[num_poly_filter]
    base_lists = base_lists[num_poly_filter]

    first_polymorphism_len = base_lists.agg(lambda l: len(l[0]))
    second_polymorphism_len = base_lists.agg(lambda l: len(l[1]))

    insersion_del_filter = (first_polymorphism_len == 1) & (second_polymorphism_len == 1)
    sys.stderr.write("{:d} insersion/del filtered.\n".format((~insersion_del_filter).sum()))

    base_info_df = base_info_df.loc[insersion_del_filter]

    # remove duplicates
    loc_str = base_info_df["chrom"] + "_" + base_info_df["start_loc"].agg(str)
    loc_count_series = loc_str.value_counts()
    dup_loc_strs = loc_count_series.index[loc_count_series > 1]

    dup_loc_filter = np.array([l not in dup_loc_strs for l in loc_str], dtype=bool)

    sys.stderr.write("{:d} duplicated locations filtered.\n".format((~dup_loc_filter).sum()))

    base_info_df = base_info_df.loc[dup_loc_filter]

    return base_info_df.copy()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Edit Fasta.")

    set_parser(parser)

    args = parser.parse_args()

    base_info_df = pd.read_csv(args.base_info_path, 
                               sep="\t", 
                               names=["chrom", 
                                      "start_loc", 
                                      "end_loc", 
                                      "name", 
                                      "score", 
                                      "strand", 
                                      "bases", 
                                      ], 
                               )
    base_info_df = filter_base_info(base_info_df)
    if args.filtered_base_info_opath:
        base_info_df.to_csv(args.filtered_base_info_opath, 
                            index=False, 
                            sep="\t", 
                            )

    base_info_chroms = base_info_df["chrom"].unique()

    mutated_records = []
    with open(args.ref_fasta, "r") as ref_fasta:
        for record in SeqIO.parse(ref_fasta, "fasta"):
            contig_name = record.id
            record_seq = record.seq

            if contig_name in base_info_chroms:
                mutated_contig = mutate_contig(record_seq, base_info_df)
            else:
                mutated_contig = record_seq
                    
            mutated_record = SeqRecord.SeqRecord(seq=mutated_contig, 
                                                 id=contig_name, 
                                                 name=record.name, 
                                                 description=record.description, 
                                                 )

            mutated_records.append(mutated_record)

    with open(args.opath, "w") as output_fasta:
        SeqIO.write(mutated_records, output_fasta, "fasta")

                
