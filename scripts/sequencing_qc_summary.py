#!/usr/bin/env python
import re
import sys
import json
import argparse

import pandas as pd

from RGTools.utils import str2bool

class SequencingQcSummary:
    '''
    Summarizing the sequencing QC results for a sample 
    based on bam flagstat and fastp json files.
    '''
    @staticmethod
    def get_fastp_fields():
        '''
        List of fields in the fastp json file that are outputted for QC.
        '''
        return ["total_reads", 
                "passed_filter_reads", 
                "percentage_passed_filter", 
                "percentage_low_quality", 
                "percentage_too_many_N", 
                "percentage_too_short", 
                "percentage_too_long", 
                ]

    @staticmethod
    def flagstat_key_map(key):
        '''
        Takes in a key in the flagstat file and returns the corresponding key in the summary.

        >>> flagstat_key_map("total_reads (QC-passed + QC-failed)")
        "total_reads"
        >>> flagstat_key_map("properly paired (100.00% : N/A)")
        "properly_paired"
        '''
        # Few hard coded mappings
        static_key_dict = {
            "in total (QC-passed reads + QC-failed reads)": "total_reads",
            "mapped": "mapped_reads",
        }

        # regex mappings to remove percentage annotations for the keys
        if match := re.match(r"(.*) \((.*): N/A\)", key):
            key = match.group(1)

        # Replace the key with the static hard coded key if it exists
        if key in static_key_dict.keys():
            key = static_key_dict[key]

        # Replace all the special characters 
        while " " in key:
            key = key.replace(" ", "_")

        while "(" in key:
            key = key.replace("(", "")
        
        while ")" in key:
            key = key.replace(")", "")

        while ">=" in key:
            key = key.replace(">=", "_geq_")
        
        while "<=" in key:
            key = key.replace("<=", "_leq_")
        
        while "<" in key:
            key = key.replace("<", "_lt_")
        
        while ">" in key:
            key = key.replace(">", "_gt_")
        
        return key
    
    @staticmethod
    def set_general_parser_arguments(parser):
        '''
        Set the general parser arguments that are used in 
        all subparsers:

        - sample: sample name for the sample to QC.
        - opath: output path for the sequencing QC summary.
        '''
        parser.add_argument("--sample",
                            help="Sample name for the sample to QC.",
                            action="append",
                            )
        
        parser.add_argument("--opath",
                            help="Output path for the sequencing QC summary. [stdout]",
                            default="stdout",
                            type=str,
                            )
        
        parser.add_argument("--output_header",
                            help="If to output header for the summary. [True]",
                            type=str2bool,
                            default=False,
                            )

    @staticmethod
    def set_trim_qc_parser(trim_qc_parser):
        SequencingQcSummary.set_general_parser_arguments(trim_qc_parser)

        trim_qc_parser.add_argument("--fastp_json_path",
                                    help="Path to fastp json file for the sample.",
                                    action="append",
                                    type=str,
                                    )
        
    @staticmethod
    def set_alignment_qc_parser(alignment_qc_parser):
        SequencingQcSummary.set_general_parser_arguments(alignment_qc_parser)

        alignment_qc_parser.add_argument("--STAR_log_path",
                                        help="Path to STAR log for the sample.",
                                        action="append",
                                        type=str,
                                        )

    @staticmethod
    def set_flagstat_qc_parser(flagstat_qc_parser):
        SequencingQcSummary.set_general_parser_arguments(flagstat_qc_parser)

        flagstat_qc_parser.add_argument("--flagstat_path",
                                        help="Path to flagstat file for the sample.",
                                        action="append",
                                        type=str,
                                        )

        flagstat_qc_parser.add_argument("--flagstat_key",
                                        help="Key in the flagstat file to extract the QC info."
                                             "(eg. total_reads, mapped_reads, properly paired, etc.)",
                                        type=str,
                                        )
        
        flagstat_qc_parser.add_argument("--output_field_name",
                                        help="Output field name for the flagstat QC info.",
                                        type=str,
                                        )
                                    
    @staticmethod
    def set_parser(parser):
        subparsers = parser.add_subparsers(dest="command")

        trim_qc_parser = subparsers.add_parser("trim_qc", help="QC the trimming step.")
        SequencingQcSummary.set_trim_qc_parser(trim_qc_parser)

        alignment_qc_parser = subparsers.add_parser("alignment_qc", help="QC the alignment step.")
        SequencingQcSummary.set_alignment_qc_parser(alignment_qc_parser)

        flagstat_qc_parser = subparsers.add_parser("flagstat_qc", help="QC flagstat files.")
        SequencingQcSummary.set_flagstat_qc_parser(flagstat_qc_parser)

    @staticmethod
    def read_flagstat(flagstat_path):
        '''
        Read a flagstat file and return the info in a dictionary.
        '''
        info_dict = {}

        line_regex = re.compile(r"(\d+) \+ (\d+) (.*)")

        with open(flagstat_path, "r") as f:
            while line := f.readline():
                match = line_regex.match(line)
                info_key = SequencingQcSummary.flagstat_key_map(match.group(3))
                mapped_reads_pass_qc = int(match.group(1))

                info_dict[info_key] = mapped_reads_pass_qc

        return info_dict
    
    @staticmethod
    def read_fastp_json(fastp_json_path):
        '''
        Read a fastp json file and return the info in a dictionary.

        Reading the following fields:
        summary: before_filtering:
        - total_reads

        filtering result:
        - passed_filter_reads
        - low_quality_reads
        - too_many_N_reads
        - too_short_reads
        - too_long_reads

        %TODO: expand the methods to read more fields
        '''
        with open(fastp_json_path, "r") as f:
            fastp_json = json.load(f)

        info_dict = {}

        info_dict["total_reads"] = fastp_json["summary"]["before_filtering"]["total_reads"]
        info_dict["passed_filter_reads"] = fastp_json["filtering_result"]["passed_filter_reads"]
        info_dict["low_quality_reads"] = fastp_json["filtering_result"]["low_quality_reads"]
        info_dict["too_many_N_reads"] = fastp_json["filtering_result"]["too_many_N_reads"]
        info_dict["too_short_reads"] = fastp_json["filtering_result"]["too_short_reads"]
        info_dict["too_long_reads"] = fastp_json["filtering_result"]["too_long_reads"]

        # statistics
        info_dict["percentage_passed_filter"] = info_dict["passed_filter_reads"] / info_dict["total_reads"]
        info_dict["percentage_low_quality"] = info_dict["low_quality_reads"] / info_dict["total_reads"]
        info_dict["percentage_too_many_N"] = info_dict["too_many_N_reads"] / info_dict["total_reads"]
        info_dict["percentage_too_short"] = info_dict["too_short_reads"] / info_dict["total_reads"]
        info_dict["percentage_too_long"] = info_dict["too_long_reads"] / info_dict["total_reads"]

        return info_dict
    
    def read_STAR_log(STAR_log_path):
        '''
        Read a STAR log file and return the info in a dictionary.
        '''
        info_dict = {}

        with open(STAR_log_path, "r") as f:
            for line in f:
                if "Number of input reads" in line:
                    info_dict["total_reads"] = int(line.split("\t")[1])
                if "Uniquely mapped reads number" in line:
                    info_dict["uniquely_mapped_reads"] = int(line.split("\t")[1])
                if "Uniquely mapped reads % " in line:
                    info_dict["percentage_uniquely_mapped_reads"] = float(line.split("\t")[1].replace("%", ""))/100
                if "Average mapped length" in line:
                    info_dict["average_mapped_length"] = float(line.split("\t")[1])
                if "Number of splices: Total" in line:
                    info_dict["total_splices"] = int(line.split("\t")[1])
                if "Number of reads mapped to multiple loci" in line:
                    info_dict["multimapper_count"] = int(line.split("\t")[1])
                if r"% of reads mapped to multiple loci" in line:
                    info_dict["percentage_multimappers"] = float(line.split("\t")[1].replace("%", ""))/100
                if "Number of reads unmapped: too many mismatches" in line:
                    info_dict["unmapped_too_many_mismatches"] = int(line.split("\t")[1])
                if r"% of reads unmapped: too many mismatches" in line:
                    info_dict["percentage_unmapped_too_many_mismatches"] = float(line.split("\t")[1].replace("%", ""))/100
                if "Number of reads unmapped: too short" in line:
                    info_dict["unmapped_too_short"] = int(line.split("\t")[1])
                if r"% of reads unmapped: too short" in line:
                    info_dict["percentage_unmapped_too_short"] = float(line.split("\t")[1].replace("%", ""))/100
                
        return info_dict
    
    @staticmethod
    def trim_qc_main(trim_qc_args):
        '''
        Main function for QC trimming results.
        '''
        result_df = pd.DataFrame(columns=["field"] + trim_qc_args.sample)

        for sample, fastp_json_path in zip(trim_qc_args.sample, trim_qc_args.fastp_json_path):
            fastp_info = SequencingQcSummary.read_fastp_json(fastp_json_path)
            for ind, field in enumerate(SequencingQcSummary.get_fastp_fields()):
                result_df.loc[ind, "field"] = "fastp_" + field
                result_df.loc[ind, sample] = fastp_info[field]

        SequencingQcSummary.output_df(result_df, 
                                      trim_qc_args.opath, 
                                      header=trim_qc_args.output_header,
                                      )
    
    @staticmethod
    def alignmenet_qc_main(alignment_qc_args):
        '''
        Main function for QC alignment results.
        '''
        result_df = pd.DataFrame(columns=["field"] + alignment_qc_args.sample)

        for sample, STAR_log_path in zip(alignment_qc_args.sample, alignment_qc_args.STAR_log_path):
            STAR_info = SequencingQcSummary.read_STAR_log(STAR_log_path)
            for ind, field in enumerate(STAR_info.keys()):
                result_df.loc[ind, "field"] = "STAR_" + field
                result_df.loc[ind, sample] = STAR_info[field]

        SequencingQcSummary.output_df(result_df,
                                      alignment_qc_args.opath,
                                      header=alignment_qc_args.output_header,
                                      )
    
    @staticmethod
    def flagstat_qc_main(flagstat_qc_args):
        '''
        Main function for flastat QC results.
        '''
        result_df = pd.DataFrame(columns=["field"] + flagstat_qc_args.sample)

        result_df.loc[0, "field"] = flagstat_qc_args.output_field_name

        for sample, flagstat_path in zip(flagstat_qc_args.sample, flagstat_qc_args.flagstat_path):
            flagstat_info = SequencingQcSummary.read_flagstat(flagstat_path)
            # Check for flagstat key existence in the flagstat file
            if not flagstat_qc_args.flagstat_key in flagstat_info.keys():
                raise ValueError(f"Flagstat key {flagstat_qc_args.flagstat_key} not found in the flagstat file {flagstat_path}".format(flagstat_qc_args.flagstat_key, 
                                                                                                                                       flagstat_path, 
                                                                                                                                       ))
            result_df.loc[0, sample] = flagstat_info[flagstat_qc_args.flagstat_key]
    
        SequencingQcSummary.output_df(result_df,
                                      flagstat_qc_args.opath,
                                      header=flagstat_qc_args.output_header,
                                      )

    @staticmethod
    def output_df(df, opath, **kwargs):
        '''
        Output the dataframe to the opath.

        Keyword arguments:
        - df: dataframe to output
        - opath: output path
        - kwargs: keyword arguments for the pd.to_csv function
        '''
        if opath == "stdout":
            df.to_csv(sys.stdout, sep="\t", index=False, **kwargs)
        else:
            df.to_csv(opath, sep="\t", index=False, **kwargs)

    @staticmethod
    def main(args) -> None:
        '''
        Main function for summarizing sequencing QC.
        
        args: 
        - sample: list of sample names
        - flagstat_path: list of paths to flagstat files
        - fastp_json_path: list of paths to fastp json files
        - opath: output path for the summary
        '''
        if args.command == "trim_qc":
            SequencingQcSummary.trim_qc_main(args)
        if args.command == "alignment_qc":
            SequencingQcSummary.alignmenet_qc_main(args)
        if args.command == "flagstat_qc":
            SequencingQcSummary.flagstat_qc_main(args)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    SequencingQcSummary.set_parser(parser)
    args = parser.parse_args()

    SequencingQcSummary.main(args)
