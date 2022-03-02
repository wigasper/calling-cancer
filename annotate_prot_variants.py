#!/usr/bin/env python
import re
import sys
import json
import logging
import logging.config
from pathlib import Path
from typing import Dict, Tuple, List

import requests

from codon_map import get_codon_map

logging.config.fileConfig("logging.conf")
logger = logging.getLogger("drivers")
logger.setLevel(logging.DEBUG)


def load_token(fp: Path = Path("oncokb_key")):
    token = None

    with open(fp, "r") as handle:
        token = handle.readline().strip()

    return token


REQ_HEADERS = {
    "application": "application/json",
    "Authorization": f"Bearer {load_token()}",
}

ALT_REQ_PREFIX = (
    "https://www.oncokb.org/api/v1/annotate/mutations/"
    "byProteinChange?referenceGenome=GRCh38&hugoSymbol="
)


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


def load_cancer_gene_list(fp: Path):
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
):
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
):
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


def known_alts_handler(known_alts_fp: Path, cancer_genes: Dict):
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
):

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


def main():
    if "-h" in sys.argv or "--help" in sys.argv:
        print(
            "Usage:\n$ python3 annotate_prot_variants.py vcfs_dir json_output_fp"
        )
        sys.exit()

    codon_map = get_codon_map()

    aa_map = get_aa_map(codon_map)

    ann_vcfs_dir = Path(sys.argv[1])

    results = {}

    for fp in ann_vcfs_dir.iterdir():
        results[fp.stem.split(".")[0]] = parse_ann_vcf(fp, aa_map)

    with open(sys.argv[2], "w") as out:
        json.dump(results, out)


if __name__ == "__main__":
    main()
