# RegGenome_script_lib
Scripts for Regulatory Genome analysis.

## Scripts

This section is a brief description of 
the usage of each scripts. For detailed 
input/output, please use the `--help` flag 
for each script.

### compute_bw_correlation.py

This program computes correlation 
between 2 bigWig files. Correlations 
are computed by aggregating signals 
in fixed genomic windows (default: 50bp).

### count_bw_sig.py

This program produce a table of signal 
counts ready for differential expression 
analysis given a list of bigWig files 
and a region file.

### gen_tss_bed_from_gtf.py

This program output a bed file with regions 
centering at TSS according to a GTF 
annotation.

### gtf_search.py

This program output a bed file with regions 
that satisfies given search criteria 
in a gtf annotation.

### merge_bw.py

This program merge the signal by base-pair 
from 2 bigWig files.

### plot_heatmap.py

This program read a count table (as produced by `merge_bw.py`) 
and plot a clustermap with seaborn.

### summarize_RSEM.py

This program converts RSEM results to 
multiple data tables (eg. TPM table, 
counts table etc.).


