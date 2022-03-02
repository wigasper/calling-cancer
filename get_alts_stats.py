#!/usr/bin/env python
import re
import sys
import json
import argparse
import logging
import logging.config
from pathlib import Path
from typing import Dict, Tuple, List

import requests

from codon_map import get_codon_map

logging.config.fileConfig("logging.conf")
logger = logging.getLogger("drivers")
logger.setLevel(logging.DEBUG)


REQ_HEADERS = {
    "application": "application/json",
    "Authorization": f"Bearer {None}",
}

ALT_REQ_PREFIX = (
    "https://www.oncokb.org/api/v1/annotate/mutations/"
    "byProteinChange?referenceGenome=GRCh38&hugoSymbol="
)


def load_token(fp: Path = Path("oncokb_key")) -> str:
    token = None

    with open(fp, "r") as handle:
        token = handle.readline().strip()

    return token


def parse_ann_alt(alteration: str, aa_map: Dict) -> Tuple[str, int, str]:
    out = None

    alteration = alteration.split(".")[1]

    match = re.search("(^\D+)(\d+)(\D.*)", alteration)

    if match:
        ref = match.group(1)
        position = int(match.group(2))
        alt = match.group(3)

        if ref in aa_map:
            ref = aa_map[ref]
        else:
            ref = "?"

        if not alt == "fs" and not alt == "*" and alt in aa_map:
            alt = aa_map[alt]

        out = (ref, position, alt)

    return out


def parse_cancer_gene_list_json(res: List[Dict]) -> Dict:
    out = {}

    for it in res:
        out[it["entrezGeneId"]] = {
            "hugo_symbol": it["hugoSymbol"],
            "refseq": it["grch38RefSeq"],
        }

    return out


def write_cancer_gene_list(cancer_genes: Dict, fp: Path):
    with open(fp, "w") as out:
        out.write("entrez_geneID\thugo_symbol\trefseq\n")
        for it, alt_ids in cancer_genes.items():
            out.write(f"{it}\t{alt_ids['hugo_symbol']}\t{alt_ids['refseq']}\n")


def load_cancer_gene_list(fp: Path) -> Dict:
    out = {}

    with open(fp, "r") as handle:
        handle.readline()
        for line in handle:
            line = line.split()
            out[line[0]] = {"hugo_symbol": line[1]}  # , "refseq": line[2]}

            if len(line) == 3:
                out[line[0]]["refseq"] = line[2]

    return out


def cancer_gene_list_handler(
    cancer_gene_list_fp: Path = Path(".cancer_gene_list"), token=load_token()
) -> Dict:
    cancer_genes = None

    if not cancer_gene_list_fp.exists():
        req_url = "https://www.oncokb.org/api/v1/utils/cancerGeneList"

        res = requests.get(req_url, headers=REQ_HEADERS)

        if not res.status_code == 200:
            raise ValueError("Cancer gene list does not exist and API call failed")

        cancer_genes = parse_cancer_gene_list_json(res.json())

        write_cancer_gene_list(cancer_genes, cancer_gene_list_fp)
    else:
        cancer_genes = load_cancer_gene_list(cancer_gene_list_fp)

    return cancer_genes


def alteration_handler(
    gene: str,
    ref: str,
    position: int,
    alt: str,
    known_alts: Dict,
    results_storage_dir: Path = Path("requests"),
) -> str:
    if not results_storage_dir.exists():
        results_storage_dir.mkdir()

    out = None

    position = str(position)

    if position in known_alts[gene] and alt in known_alts[gene][position]:
        out = known_alts[gene][position][alt]
    else:
        req_url = f"{ALT_REQ_PREFIX}{gene}&alteration={ref}{position}{alt}"

        res = requests.get(req_url, headers=REQ_HEADERS)

        if res.status_code == 200:
            this_res = res.json()
            with open(
                results_storage_dir.joinpath(Path(f"{gene}|{ref}{position}{alt}")), "w"
            ) as out:
                json.dump(this_res, out)

            if position not in known_alts[gene]:
                known_alts[gene][position] = {}

            known_alts[gene][position][alt] = this_res["oncogenic"]
            out = this_res["oncogenic"]

        else:
            logger.info(f"API call failure for {gene} {ref}{position}{alt}")

    return out


def known_alts_handler(known_alts_fp: Path, cancer_genes: Dict) -> Dict:
    out = None

    if known_alts_fp.exists():
        with open(known_alts_fp, "r") as handle:
            out = json.load(handle)
    else:
        out = {}
        for _entrez_uid, alts in cancer_genes.items():
            out[alts["hugo_symbol"]] = {}

    return out


def parse_ann_vcf(
    fp: Path,
    aa_map: Dict,
    cancer_gene_list_fp: Path = Path(".cancer_gene_list"),
    known_alts_fp: Path = Path(".memoized_alts"),
) -> Dict:

    cancer_genes = cancer_gene_list_handler(cancer_gene_list_fp)
    known_alts = known_alts_handler(known_alts_fp, cancer_genes)

    cancer_gene_symbols = set(
        [vals["hugo_symbol"] for _k, vals in cancer_genes.items()]
    )

    onco_stats = {}

    variants_of_interest = [
        "missense_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
    ]

    with open(fp, "r") as handle:
        for line in handle:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                info = line[7]

                info = info.split(";")

                info = [it for it in info if it.startswith("ANN")]

                for it in info:
                    it = it.split(",")

                    it = it[0]
                    it = it.split("|")

                    variant_types = it[1].split("&")
                    gene = it[3]

                    interesting = False

                    if gene in cancer_gene_symbols:
                        # print("\t".join(line))
                        for var_type in variant_types:
                            if var_type in variants_of_interest:
                                interesting = True

                    if interesting:
                        alteration = it[10]

                        parse_res = parse_ann_alt(alteration, aa_map)

                        if parse_res:
                            ref, position, alt = parse_res
                            oncogenic_status = alteration_handler(
                                gene, ref, position, alt, known_alts
                            )

                            if oncogenic_status:
                                if oncogenic_status not in onco_stats:
                                    onco_stats[oncogenic_status] = []

                                onco_stats[oncogenic_status].append(
                                    f"{gene}|{ref}{position}{alt}"
                                )

    with open(known_alts_fp, "w") as out:
        json.dump(known_alts, out)

    return onco_stats


# map 3 letter AAs to 1 letter AAs
def get_aa_map(codon_map: Dict) -> Dict:
    aa_map = {}

    for it in codon_map:
        aa_map[codon_map[it][3]] = codon_map[it][1]

    aa_map["Ter"] = "*"

    return aa_map


def load_normal_counts(fp: Path) -> List[int]:
    counts = []

    with open(fp, "r") as handle:
        for line in handle:
            counts.append(int(line.strip()))

    return sorted(counts)


def get_percentile_stats(counts: Dict, normal_counts: List[int]) -> Dict[str, float]:
    out = {it: None for it in counts}

    for sample_uid, count in counts.items():
        pctile = len([it for it in normal_counts if it < count]) / len(normal_counts)

        out[sample_uid] = pctile

    return out


def get_stats(results: Dict, normal_counts_fp: Path) -> Dict[str, float]:
    annotations_of_interest = ["Predicted Oncogenic", "Likely Oncogenic", "Oncogenic"]

    normal_counts = load_normal_counts(normal_counts_fp)

    counts_per_cell = {}

    for sample_uid in results:
        counts_per_cell[sample_uid] = 0
        these_annotations = [
            a for a in results[sample_uid] if a in annotations_of_interest
        ]

        for annotation in these_annotations:
            for alt in results[sample_uid][annotation]:
                if not alt.startswith("HLA"):
                    counts_per_cell[sample_uid] += 1

    stats = get_percentile_stats(counts_per_cell, normal_counts)

    return stats


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        help="Directory containing SnpEff-annotated VCFs to be processed",
        required=True,
    )
    parser.add_argument(
        "-o", "--output", help="Output file name", default="driver_stats.tsv"
    )
    parser.add_argument(
        "-k",
        "--key",
        help="OncoKB API key: a file with a single line containing the key",
        default="oncokb_key",
    )
    parser.add_argument(
        "-n", "--normal-counts", help="Path to normal counts file", required=True
    )

    args = parser.parse_args()

    ann_vcfs_dir = Path(args.input)
    out_fp = Path(args.output)
    oncokb_key_fp = Path(args.key)
    normal_counts_fp = Path(args.normal_counts)

    global REQ_HEADERS
    REQ_HEADERS = {
        "application": "application/json",
        "Authorization": f"Bearer {load_token(oncokb_key_fp)}",
    }

    codon_map = get_codon_map()

    aa_map = get_aa_map(codon_map)

    results = {}

    for fp in ann_vcfs_dir.iterdir():
        results[fp.stem.split(".")[0]] = parse_ann_vcf(fp, aa_map)

    stats = get_stats(results, normal_counts_fp)

    with open(out_fp, "w") as out:
        for sample_uid, pctile in stats.items():
            out.write(f"{sample_uid}\t{pctile}\n")


if __name__ == "__main__":
    main()
