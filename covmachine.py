#!/usr/bin/env python3
import re
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, Tuple

import gtf_parser


def parse_alt(alteration: str) -> Tuple[str, int, str]:
    out = None

    # remove "driver"
    alteration = alteration.split()[0]

    match = re.search("(^\D)(\d+)(\D.*)", alteration)

    if match:
        ref = match.group(1)
        position = int(match.group(2))
        alt = match.group(3)

        if "fs" in alt:
            alt = "fs"

        out = (ref, position, alt)

    return out


def get_all_alterations(alterations_fp: Path) -> Dict:
    out = {}

    alterations = None

    with open(alterations_fp, "r") as handle:
        alterations = json.load(handle)

    for sample_uid, alts in alterations.items():
        if alts:
            # classif is an annotation category: 'Likely Oncogenic',
            # 'Inconclusive', etc
            for classif in alterations[sample_uid]:
                for alt in alterations[sample_uid][classif]:
                    gene = alt.split("|")[0]

                    if gene not in out:
                        out[gene] = {}

                    parse_res = parse_alt(alt.split("|")[1])

                    if parse_res:
                        ref, position, alt = parse_res
                        alter_dict = {"ref": ref, "alt": alt}

                        if (
                            position in out[gene]
                            and alter_dict not in out[gene][position]
                        ):
                            out[gene][position].append(alter_dict)
                        elif position not in out[gene]:
                            out[gene][position] = [alter_dict]

    return out


# Get the genomic positions corresponding to driver residues
# First, get the important genes and resides from the alterations_fp
# Then, parse the GTF to find the genomic positions for all these residues
def get_driver_positions(alterations_fp: Path, gtf_fp: Path) -> Dict:
    driver_positions = {}

    alterations = get_all_alterations(alterations_fp)
    driver_positions = {gene: {} for gene in alterations}

    for gene in alterations:
        for residue in alterations[gene]:
            driver_positions[gene][residue] = []

    features = gtf_parser.parse_gtf(gtf_fp, set(driver_positions.keys()))

    position_lookups = {}

    for feat in features:
        gtf_parser.process_transcript(feat, features, position_lookups)

    for gene in alterations:
        for residue in alterations[gene]:
            try:
                positions = position_lookups[gene][residue]
            except KeyError as e:
                # TODO: TP53 residue 394, make some logic to deal
                print(repr(e))
                print(f"KeyError looking for position_lookups[{gene}][{residue}]")
                print("Check validity")
            driver_positions[gene][residue].extend(positions)

    return driver_positions


# Gets coverage for the desired driver positions by reading sorted
# genotype likelihoods (mpileup output) through stdin
def get_coverages(driver_positions: Dict) -> Dict:
    coverages = {gene: {} for gene in driver_positions}

    for gene in driver_positions:
        for residue in driver_positions[gene]:
            coverages[gene][residue] = {}

    for gene in driver_positions:
        for residue in driver_positions[gene]:
            for chromosome, base in driver_positions[gene][residue]:
                coverages[gene][residue][base] = 0

    locs_of_interest = {}

    for gene in driver_positions:
        for residue in driver_positions[gene]:
            for chromosome, base_pos in driver_positions[gene][residue]:
                locs_of_interest[f"{chromosome}|{base_pos}"] = (gene, residue)

    print("Pulling from stdin")
    for line in sys.stdin:
        if not line.startswith("#"):
            line = line.split()
            chromosome = line[0]
            position = int(line[1])

            if f"{chromosome}|{position}" in locs_of_interest:
                gene, residue = locs_of_interest[f"{chromosome}|{position}"]

                if line[7].split(";")[0] != "INDEL":
                    depth = int(line[7].split(";")[0].split("=")[1])

                    coverages[gene][residue][position] = depth

    print("Got EOF, all done")

    return coverages


# Writes output
def write_coverages(coverages: Dict, out_fp: Path):
    with open(out_fp, "w") as out:
        for gene in coverages:
            for residue in coverages[gene]:
                for base in coverages[gene][residue]:
                    depth = coverages[gene][residue][base]
                    out.write(f"{gene}\t{residue}\t{base}\t{depth}\n")


def main():
    module = sys.argv[1]

    parser = argparse.ArgumentParser()
    # parser.add_argument("-a", "--alterations-fp", default="alterations_across_samples.tsv")
    parser.add_argument(
        "-a",
        "--alterations-fp",
        help="Path to JSON detailing alterations for all samples",
    )
    parser.add_argument("-g", "--gtf", default="hg38.ncbiRefSeq.gtf")
    parser.add_argument("-d", "--drivers", default="driver_positions.json")

    driver_positions = None

    if module == "init":
        args = parser.parse_args(sys.argv[2:])

        driver_positions = get_driver_positions(
            Path(args.alterations_fp), Path(args.gtf)
        )

        with open(Path(args.drivers), "w") as out:
            json.dump(driver_positions, out)

    if module == "analyze":
        parser.add_argument("-o", "--out", required=True)

        args = parser.parse_args(sys.argv[2:])

        if Path(args.drivers).exists():
            with open(Path(args.drivers), "r") as handle:
                driver_positions = json.load(handle)

        else:
            print("Driver positions file not found, figuring it out")
            driver_positions = get_driver_positions(
                Path(args.alterations_fp), Path(args.gtf)
            )
            print("Got driver positions")

        coverages = get_coverages(driver_positions)

        print(f"Writing output at {args.out}")
        write_coverages(coverages, Path(args.out))


if __name__ == "__main__":
    main()
