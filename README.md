# RegGenome_script_lib
Scripts for Regulatory Genome analysis.

## Testing

Before running tests, sample data 
of large size should be placed under 
`large_sample_data/`.

- `hg38.fa`: UCSC hg38 reference genome.

Test scripts are under `tests/`. 
Use the following command to 
run test for each scripts.

```{bash}
SCRIPT=#test script to run
python -m unittest ${SCRIPT}
```

## Scripts

For detailed input/output of each script, 
please use the `--help` flag for each script.

## Key concepts

Here are a few concepts the scripts center on.

### Exogeneous sequences

A set of sequences that are independent from the 
reference genome. In our convention, exogeneous 
sequences are stored as a combination of a few 
files.

- `fasta`: A fasta file stores the sequence information. 
           The boundry of the sequence do not necessarily 
           represents the boundry of exogeneous elements, 
           as they are defined by the bed region.
- `bed`: A bed file stores the region information, aka the 
         element boundry.
- `npy`: Zero or more npy file that store the necessary information, 
         storing numpy array of data in the same order 
         of the bed file.

As constricted by the numpy array data structure, informations 
are of the same length. For example, when storing data such as 
singal track (e.g. PROcap signals), elements will be padded 
to the same length. 

Although we usually put one fasta entry for each element, 
elements are defined by the bed file and can be defined on 
the same fasta entry. 

### Genomic elements

A special type of `exogeneous sequences`. The 
`fasta` file for genomic regions are restricted to 
a reference genome.

### Mutagenesis Results

Another special type of `exogeneous sequences` and can 
be processed by `exogeneous_tool`. This data type has an 
additional constraint that each bed entry must correspond 
to a fasta entry, so that results can be processed 
using pattern matching with regular expressions.
