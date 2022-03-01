#!/usr/bin/env python3

# Just returns the codon map
def get_codon_map():
    codon_map = {
        "att": {1: "I", 3: "Ile"},
        "atc": {1: "I", 3: "Ile"},
        "ata": {1: "I", 3: "Ile"},
        "ctt": {1: "L", 3: "Leu"},
        "ctc": {1: "L", 3: "Leu"},
        "cta": {1: "L", 3: "Leu"},
        "ctg": {1: "L", 3: "Leu"},
        "tta": {1: "L", 3: "Leu"},
        "ttg": {1: "L", 3: "Leu"},
        "gtt": {1: "V", 3: "Val"},
        "gtc": {1: "V", 3: "Val"},
        "gta": {1: "V", 3: "Val"},
        "gtg": {1: "V", 3: "Val"},
        "ttt": {1: "F", 3: "Phe"},
        "ttc": {1: "F", 3: "Phe"},
        "atg": {1: "M", 3: "Met"},
        "tgt": {1: "C", 3: "Cys"},
        "tgc": {1: "C", 3: "Cys"},
        "gct": {1: "A", 3: "Ala"},
        "gcc": {1: "A", 3: "Ala"},
        "gca": {1: "A", 3: "Ala"},
        "gcg": {1: "A", 3: "Ala"},
        "ggt": {1: "G", 3: "Gly"},
        "ggc": {1: "G", 3: "Gly"},
        "gga": {1: "G", 3: "Gly"},
        "ggg": {1: "G", 3: "Gly"},
        "cct": {1: "P", 3: "Pro"},
        "ccc": {1: "P", 3: "Pro"},
        "cca": {1: "P", 3: "Pro"},
        "ccg": {1: "P", 3: "Pro"},
        "act": {1: "T", 3: "Thr"},
        "acc": {1: "T", 3: "Thr"},
        "aca": {1: "T", 3: "Thr"},
        "acg": {1: "T", 3: "Thr"},
        "tct": {1: "S", 3: "Ser"},
        "tcc": {1: "S", 3: "Ser"},
        "tca": {1: "S", 3: "Ser"},
        "tcg": {1: "S", 3: "Ser"},
        "agt": {1: "S", 3: "Ser"},
        "agc": {1: "S", 3: "Ser"},
        "tat": {1: "Y", 3: "Tyr"},
        "tac": {1: "Y", 3: "Tyr"},
        "tgg": {1: "W", 3: "Trp"},
        "caa": {1: "Q", 3: "Gln"},
        "cag": {1: "Q", 3: "Gln"},
        "aat": {1: "N", 3: "Asn"},
        "aac": {1: "N", 3: "Asn"},
        "cat": {1: "H", 3: "His"},
        "cac": {1: "H", 3: "His"},
        "gaa": {1: "E", 3: "Glu"},
        "gag": {1: "E", 3: "Glu"},
        "gat": {1: "D", 3: "Asp"},
        "gac": {1: "D", 3: "Asp"},
        "aaa": {1: "K", 3: "Lys"},
        "aag": {1: "K", 3: "Lys"},
        "cgt": {1: "R", 3: "Arg"},
        "cgc": {1: "R", 3: "Arg"},
        "cga": {1: "R", 3: "Arg"},
        "cgg": {1: "R", 3: "Arg"},
        "aga": {1: "R", 3: "Arg"},
        "agg": {1: "R", 3: "Arg"},
        "taa": {1: "*", 3: "***"},
        "tag": {1: "*", 3: "***"},
        "tga": {1: "*", 3: "***"},
    }

    return codon_map
