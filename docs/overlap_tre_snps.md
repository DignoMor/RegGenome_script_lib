# Overlap TRE SNPs

The script overlaps TRE with risk SNPs 
and report a list of candidate elements 
for reporter assay.

## Input

### TRE and SNPs

TRE and SNPs are inputted in the 
format of bed files through 
the `--tre_inpath` and 
`--snp_inpath` flag.

### filter options

There are a few filter options so that 
the output is compatible for element 
synthesis.

* `max_length`: The maximum length 
for the final elements.

### Output options

The output path is defined by the 
`--opath` flags. A few output files 
will be produced:

* `${project_name}.tre.bed`: This files includes a list of TREs that are ready for reporter assay.
* `${project_name}.tre.filtered.bed`: This file includes a list of TREs that are filtered. Iinformation regarding the reason for filtering can be retrieved from `${project_name}.info.json`.
* `${project_name}.id_info.tsv`: This provide a table mapping TRE ids to a meaningful element name.
* `${project_name}.info.json`: This file contains additional infos regarding the overlapping.
