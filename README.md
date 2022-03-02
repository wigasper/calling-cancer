# README

The primary results from the paper can be reproduced with `pipe.sh`

We have also included a script, `get_alts_stats.py`, that can be used to generate
percentile statistics for cells, describing where they lie on the normal cell
putative driver alteration distributions. 

`get_alts_stats.py` takes a directory, containing VCF files as input, an OncoKB
API key file, a file containing normal cell counts (given by `normal_counts_vec`), 
and optionally an output file path:

```bash
usage: get_alts_stats.py [-h] -i INPUT [-o OUTPUT] [-k KEY] -n NORMAL_COUNTS

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Directory containing SnpEff-annotated VCFs to be processed
  -o OUTPUT, --output OUTPUT
                        Output file name
  -k KEY, --key KEY     OncoKB API key: a file with a single line containing the key
  -n NORMAL_COUNTS, --normal-counts NORMAL_COUNTS
                        Path to normal counts file
```

The OncoKB API key should be present in a single line file just containing the key:

```bash
$ cat oncokb_key 
abcdefgh-1234-4321
```

The output file is a tab-delimited file with two columns, one for sample UID (the
stem of the corresponding VCF), and one for the percentile statistic.

The normal counts vector was generated according to the methodology used in the 
paper, by randomly selecting cells from `GSE130473`, `GSE81547`, `GSE109822`,
and `SRP271375`.

## File descriptions
* `annotate_prot_variants.py` - Annotates SnpEff output VCFs with OncoKB data
* `codon_map.py` - Just the codon map, separate file for readability of the others
* `covmachine.py` - Pipe genome pileups into this to track coverage at bases of interest
* `get_alts_stats.py` - Generate statistics for cells from SnpEff-annotated VCFs
* `gtf_parser.py` - Parses a GTF to extract bases corresponding to residues of interest
* `pipe.sh` - pipeline to reproduce TNBC dataset results
* `samples.tsv` - TNBC dataset SRA UIDs and sample UIDs

## OncoKB API key requirement

Reproducing this work or running the annotation script will requires an 
OncoKB API key, which can be obtained at [OncoKB](https://www.oncokb.org/apiAccess).

`annotate_prot_variants.py` will by look for this in a single-line text
file located at `./oncokb_key`. 

