
import re
import sys
import json
import argparse

import pandas as pd

class SequencingQcSummary:
    '''
    Summarizing the sequencing QC results for a sample 
    based on bam flagstat and fastp json files.
    '''
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
    def set_parser(parser):
        parser.add_argument("--sample", 
                            help="Sample name for the sample to QC.", 
                            action="append",
                            )
        
        parser.add_argument("--flagstat_path",
                            help="Path to flagstat text file for the sample.",
                            action="append",
                            )
        
        parser.add_argument("--fastp_json_path",
                            help="Path to fastp json file for the sample.",
                            action="append",
                            )
        
        parser.add_argument("--flagstat_fields", 
                            help="Fields to extract from the flagstat file. [total_reads]",
                            action="append",
                            default=["total_reads"],
                            )
        
        parser.add_argument("--fastp_fields",
                            help="Fields to extract from the fastp json file. [total_reads]",
                            action="append",
                            default=["total_reads"],
                            )

        parser.add_argument("--opath",
                            help="Output path for the sequencing QC summary. [stdout]",
                            default="stdout",
                            type=str,
                            )

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

    @staticmethod
    def main(args) -> None:
        '''
        Main function ofr summarizing sequencing QC.
        
        args: 
        - sample: list of sample names
        - flagstat_path: list of paths to flagstat files
        - fastp_json_path: list of paths to fastp json files
        - opath: output path for the summary
        '''
        result_df = pd.DataFrame(columns=["sample"] + \
                                 ["flagstat_" + f for f in args.flagstat_fields] + \
                                 ["fastp_" + f for f in args.fastp_fields])
        result_df.set_index("sample", inplace=True)

        for sample, flagstat_path, fastp_json_path in zip(args.sample, args.flagstat_path, args.fastp_json_path):
            flagstat_info = SequencingQcSummary.read_flagstat(flagstat_path)
            fastp_info = SequencingQcSummary.read_fastp_json(fastp_json_path)

            for field in args.flagstat_fields:
                result_df.loc[sample, "flagstat_" + field] = flagstat_info[field]
            
            for field in args.fastp_fields:
                result_df.loc[sample, "fastp_" + field] = fastp_info[field]
        
        result_df.reset_index(inplace=True)

        if args.opath == "stdout":
            result_df.to_csv(sys.stdout, sep="\t", index=False)
        else:
            result_df.to_csv(args.opath, sep="\t", index=False)


    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    SequencingQcSummary.set_parser(parser)
    args = parser.parse_args()

    SequencingQcSummary.main(args)
