# Updating in progress

## File descriptions
* `pipe.sh` - pipeline to reproduce TNBC dataset results
* `samples.tsv` - TNBC dataset SRA UIDs and sample UIDs
* `annotate_prot_variants.py` - Annotates SnpEff output VCFs with OncoKB data
* `codon_map.py` - Just the codon map, separate file for readability of the others
* `covmachine.py` - Pipe genome pileups into this to track coverage at bases of interest
* `gtf_parser.py` - Parses a GTF to extract bases corresponding to residues of interest

## OncoKB API key requirement

Reproducing this work or running the annotation script will requires an 
OncoKB API key, which can be obtained at [OncoKB](https://www.oncokb.org/apiAccess).

`annotate_prot_variants.py` will by default look for this in a single-line text
file located at `./oncokb_key`, but the path can also be specified as the last
argument in the `annotate_prot_variants.py` call 
(see `$python3 annotate_prot_variants.py -h`). If using an alternate key file
path and running the `pipe.sh` to reproduce this work, you will need to specify
this at line 164 of `pipe.sh` by adding the key file path as the last argument.


