#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Script to run IgBLAST on a given input FASTA file and output the results in TSV format.
import os
import subprocess
import argparse
import sys
import pandas as pd

def run_igblast(input_fasta, output_tsv, reference_dir):
    """
    Run IgBLAST on the input FASTA file and save the results to the output TSV file.

    Parameters:
    - input_fasta: Path to the input FASTA file containing sequences to be analyzed.
    - output_tsv: Path to the output TSV file where results will be saved.
    - reference_db: Path to the IgBLAST reference database.
    """
    command_igblastn = [
        "igblastn",
        "-germline_db_V",
        f"{reference_dir}/database/imgt_human_ig_v",
        "-germline_db_D",
        f"{reference_dir}/database/imgt_human_ig_d",
        "-germline_db_J",
        f"{reference_dir}/database/imgt_human_ig_j",
        "-query",
        input_fasta,
        "-organism",
        "human",
        "-auxiliary_data",
        f"{reference_dir}/optional_file/human_gl.aux",
        "-show_translation",
        "-outfmt",
        "19",
        "-out",
        output_tsv,
    ]
    try:
        subprocess.run(command_igblastn, check=True)
        print(f"IgBLAST completed successfully. Results saved to {output_tsv}")
    except subprocess.CalledProcessError as e:
        print(f"Error running IgBLAST: {e}", file=sys.stderr)
        sys.exit(1)

def convert_tsv_to_fasta(tsv_file, fasta_file):
    """
    Convert IgBLAST TSV output to FASTA format.

    Parameters:
    - tsv_file: Path to the IgBLAST TSV output file.
    - fasta_file: Path to the output FASTA file.
    """
    tsv = pd.read_csv(tsv_file, sep='\t', comment='#')
    with open(fasta_file, 'w') as fasta:
        for index, row in tsv.iterrows():
            fasta.write(f">{row['sequence_id']}\n{row['sequence']}\n")

    print(f"Converted {tsv_file} to {fasta_file}")

def keep_metadata(original_tsv, igblast_tsv, output_tsv):
    """
    Merge metadata from the original TSV file with IgBLAST results.

    Parameters:
    - original_tsv: Path to the original TSV file with metadata.
    - igblast_tsv: Path to the IgBLAST TSV output file.
    - output_tsv: Path to the final output TSV file with merged data.
    """
    original_df = pd.read_csv(original_tsv, sep='\t')
    igblast_df = pd.read_csv(igblast_tsv, sep='\t', comment='#')

    # check column names that are not present in igbast_df
    columns_to_merge = [col for col in original_df.columns if col not in igblast_df.columns]
    merged_df = pd.merge(igblast_df, original_df[columns_to_merge + ['sequence_id']], on='sequence_id', how='right')
    merged_df.to_csv(output_tsv, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Run IgBLAST on input TSV file.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("-r", "--reference", required=True, help="Path to IgBLAST reference database directory")

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Input file {args.input} does not exist.", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.reference):
        print(f"Reference directory {args.reference} does not exist.", file=sys.stderr)
        sys.exit(1)

    print(f"Running IgBLAST on {args.input} with reference {args.reference}")
    # Convert tsv to fasta
    fasta_file = args.input.replace(".tsv", ".fasta")
    convert_tsv_to_fasta(args.input, fasta_file)

    # Run IgBLAST
    igblast_output = args.output.replace(".tsv", "_igblastintermediate.tsv")
    run_igblast(fasta_file, igblast_output, args.reference)

    # Keep metadata in igblast output
    keep_metadata(args.input, igblast_output, args.output)

    print(f"Final IgBLAST results with metadata saved to {args.output}")

if __name__ == "__main__":
    main()

