#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Preprocess repertoire file by subsampling and translating sequences.
This script reads a repertoire TSV file, drops rows with incomplete CDR1 sequences,
generates a new nucleotide sequence column without FWR1, translates the sequences to amino acids. The output is saved to a new TSV file.

Usage:
python preprocess_repertoire.py --input_repertoire <input_file.tsv> --outname <output_file.tsv>
"""

import pandas as pd
from rpy2.robjects import pandas2ri, default_converter
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import r
import argparse
import os

alakazam = importr('alakazam')

# parameters
parser = argparse.ArgumentParser(description="Preprocess repertoire file")
parser.add_argument("--input_repertoire", type=str, required=True, help="Path to input repertoire TSV file")
# Do not subsample here as we want to use all sequences for the selection step
#parser.add_argument("--subsample_size", type=int, default=100000, help="Number of rows to subsample")
parser.add_argument("--outname", type=str, default=None, help="Output file name (optional)")
args = parser.parse_args()

input_repertoire = args.input_repertoire
#subsample_size = args.subsample_size
if args.outname:
    outname = args.outname
else:
    outname = os.path.basename(input_repertoire).replace('.tsv', '_subsampled.tsv')


repertoire = pd.read_csv(input_repertoire, sep="\t", header=0)


repertoire.columns


if "sequence_vdj_aa" in repertoire.columns:
    print("Existing column 'sequence_vdj_aa' will be dropped.")
    repertoire = repertoire.drop(columns=["sequence_vdj_aa"])


print("Repertoire size:", repertoire.shape)
original_size = repertoire.shape[0]

# drop rows where "cdr1" column starts with "."
repertoire = repertoire[repertoire["cdr1"].notnull()]
repertoire = repertoire[~repertoire["cdr1"].str.startswith(".")]

# drop rows with incomplete fwr4 sequences
repertoire = repertoire[repertoire["fwr4"].notnull()]
repertoire = repertoire[~repertoire["fwr4"].str.endswith("-")]

rows_dropped = original_size - repertoire.shape[0]

print("Dropped sequences with incompleted cdr1 or fwr4:", rows_dropped)


# generate new sequence_nt column by concatenating cdr1, fwr2, cdr2, fwr3, cdr3, fwr4
repertoire["sequence_nt_nofwr1"] = (
    repertoire["cdr1"]
    + repertoire["fwr2"]
    + repertoire["cdr2"]
    + repertoire["fwr3"]
    + repertoire["cdr3"]
    + repertoire["fwr4"]
)
repertoire["sequence_nt_nofwr1"] = repertoire["sequence_nt_nofwr1"].str.replace("-", "").str.replace(".", "")


translation = list(alakazam.translateDNA(repertoire["sequence_nt_nofwr1"].tolist()))


repertoire["sequence_vdj_aa"] = translation


#subset = repertoire.sample(n=subsample_size, random_state=42)


repertoire.to_csv(outname, sep="\t", index=False)


