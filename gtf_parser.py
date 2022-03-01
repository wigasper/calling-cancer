#!/usr/bin/env python3
import re
import time
import sys
from pathlib import Path
from typing import Set, Tuple, List, Dict


# Processes a transcript using the features extracted from the GTF
# Creates a lookup table for genomic positions corresponding to 
# every residue in the coding sequence
# Modifies position_lookups in place (borrow here)
def process_transcript(feature: str, features: Dict, position_lookups: Dict):
    gene = features[feature]["gene_id"]
    chromosome = features[feature]["chromosome"]
    try:
        strand = features[feature]["transcript"][0][2]
    except:
        print(feature)
        print(features[feature])
        sys.exit()

    for cds in features[feature]["cds"]:
        if cds[2] != strand:
            raise ValueError(f"Strand disagreement {feature}")

    if gene not in position_lookups:
        position_lookups[gene] = {}

    # These variables make sure that we are processing coding
    # seqs appropriately for the strand that they are located
    # on
    step_modifier = 1
    cds_order_mod = (0, 1)
    
    # creative mathematical way to map 1 to (0,1) and -1 to (1, 0) ? 
    if strand == "-":
        step_modifier = -1
        cds_order_mod = (1, 0)

    idx = 0
    for cds in features[feature]["cds"][::step_modifier]:
        for base_pos in range(cds[cds_order_mod[0]], cds[cds_order_mod[1]] + step_modifier, step_modifier):
            residue_posit = (idx // 3) + 1
            if (
                residue_posit in position_lookups[gene]
                and (chromosome, base_pos)
                not in position_lookups[gene][residue_posit]
            ):
                position_lookups[gene][residue_posit].append((chromosome, base_pos))
            else:
                position_lookups[gene][residue_posit] = [(chromosome, base_pos)]
            idx += 1


# Parses a GTF file
# Runs through every line, when a gene of interest is found
# it records unique identifiers, locations, and CDS data
def parse_gtf(gtf_fp: Path, features_of_interest: Set) -> Dict:
    feature_data = {}
    
    with open(gtf_fp, "r") as handle:
        for line in handle:
            line = line.split("\t")
            info = line[8].split(";")
            gene_id = info[0].split()[-1].strip('"')
            transcript_id = info[1].split()[-1].strip('"')
            
            #forward_condition = (gene_id in features_of_interest and feature_type == "genes") #or 
                #(transcript_id in features_of_interest and feature_type == "transcripts"))
            
            if gene_id in features_of_interest or transcript_id in features_of_interest:#forward_condition:
                chromosome = line[0]

                if line[2] in ["CDS", "transcript", "3UTR"]:
                    feature_type = line[2].lower()
                    start = int(line[3])
                    end = int(line[4])
                    strand = line[6]

                    if feature_type == "cds":
                        frame = int(line[7])
                    else:
                        frame = None
                    
                    transcript_id = info[1].split()[-1].strip('"')
                    
                    if transcript_id in feature_data:
                        feature_data[transcript_id][feature_type].append((start, end, strand, frame))
                    else:
                        feature_data[transcript_id] = {
                                    "transcript": [],
                                    "3utr": [],
                                    "cds": [],
                                    "gene_id": gene_id,
                                    "chromosome": chromosome
                                }
                        feature_data[transcript_id][feature_type].append((start, end, strand, frame))
    
    return feature_data




# Fundamentally a function solely for testing, but may be useful
# in the future. Loads genes and residues of interest from a text
# file of
# Gene Residue
def load_genes_positions(fp: Path) -> Tuple[List, List]:
    genes_of_interest = []
    positions = []

    with open(fp, "r") as handle:
        for line in handle:
            line = line.strip("\n").split()
            genes_of_interest.append(line[0])
            position = re.search("(\d+)", line[1]).group(1)
            positions.append(int(position))

    return genes_of_interest, positions


# Main routine is really only currently used for testing
# Functions are really meant to be used by covmachine.py
def main():
    genes_positions_fp = Path(sys.argv[1])
    gtf_fp = Path(sys.argv[2])
    output_fp = Path(sys.argv[3])

    genes_of_interest, positions = load_genes_positions(genes_positions_fp)
    # leaving for later, probably want to do something fun here
    #    raw_input = []

    #    for line in sys.stdin:
    #        raw_input.append(line.strip())

    #    genes_of_interest, positions = zip(*[it.split() for it in raw_input])
    # genes_of_interest = set(["KRAS", "TP53"])
    
    start = time.perf_counter()

    features = parse_gtf(gtf_fp, genes_of_interest)
    
    position_lookups = {}

    for feature in features:
        process_transcript(feature, features, position_lookups)

    for gene, position in zip(genes_of_interest, positions):
        print(f"{gene}\t{position}\t{position_lookups[gene][int(position)]}")

if __name__ == "__main__":
    main()
