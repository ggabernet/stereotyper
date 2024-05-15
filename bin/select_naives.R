#!/usr/bin/env Rscript
# Written by Gisela Gabernet and released under the MIT license (2020).


# Select naive sequences from a repertoire:
# Arguments:
#   --repertoire    Tabulated data in AIRR (TSV) format with clonal assignments and germline assignments.
#   --target        Amino acid sequence to be used as target for the convergence process.
#   --target_vgene      Target V gene.
#   --target_jgene      Target J gene.
#   --target_cdr3       Target CDR3 amino acid sequence.
#   --max_mu_freq       Maximum mutation frequency to consider for selecting naives.
#   --selection_number  Number of sequences to select. Default is 1.
#   --random_seed       Random seed. Default is 1234.
#   --outname       Output file name.
#   -h  Display help.
# Example: ./select_naives.R --repertoire AIRR_repertoire.tsv --target "QVS" --target_vgene "IGHV1-69" --target_jgene "IGHJ6" --outname "selected_naives.fasta"

# Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(seqinr))

# Define commmandline arguments
opt_list <- list(
    make_option(c("--repertoire"), default=NULL,
                help="Input repertoire .tsv file after clone definition and germline definition."),
    make_option(c("--target"), default=NULL,
                help="Amino acid sequence to be used as target for the convergence process."),
    make_option(c("--target_vgene"), default=NULL,
                help="Target V gene."),
    make_option(c("--target_jgene"), default=NULL,
                help="Target J gene."),
    make_option(c("--target_cdr3"), default=NULL,
                help="Target CDR3 amino acid sequence."),
    make_option(c("--max_mu_freq"), default = NULL, help="Maximum mutation frequency to consider for selecting naives."),
    make_option(c("--selection_number"), default = 1, help="Number of sequences to select."),
    make_option(c("--random_seed"), default = 42, help="Random seed."),
    make_option(c("--outname"), help="Output file name.")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# Set random seed
set.seed(as.numeric(opt$random_seed))

# Read repertoire
repertoire <- airr::read_rearrangement(opt$repertoire)
target <- opt$target
target_vgene <- opt$target_vgene
target_jgene <- opt$target_jgene
target_cdr3 <- opt$target_cdr3
target_cdr3_length <- nchar(target_cdr3)
max_mu_freq <- as.numeric(opt$max_mu_freq)
selection_number <- as.numeric(opt$selection_number)

# Check for needed columns
if (!("v_call" %in% colnames(repertoire)) | !("j_call" %in% colnames(repertoire)) | !("sequence_alignment" %in% colnames(repertoire)) | !("sequence" %in% colnames(repertoire)) | !("subject_id" %in% colnames(repertoire))){
    if (!("cdr3_aa" %in% colnames(repertoire)) & !("junction_aa" %in% colnames(repertoire))){
        stop("The repertoire file must contain columns 'v_call', 'j_call', 'sequence_alignment', 'sequence' and 'subject_id'.\n
        Additionally, it must contain the column 'cdr3_aa' or 'junction_aa'.\n")
    }
}

if (!("cdr3_aa" %in% colnames(repertoire))){
    repertoire$cdr3_aa <- substr(repertoire$junction_aa, 2, nchar(repertoire$junction_aa)-1)
}

# Get V gene and J gene
repertoire$v_gene <- alakazam::getGene(repertoire$v_call)
repertoire$j_gene <- alakazam::getGene(repertoire$j_call)

print(repertoire)

# Select compatible sequences
selected_naives <- repertoire %>%
    dplyr::mutate(cdr3_aa_length = nchar(cdr3_aa)) %>%
    dplyr::filter(v_gene == target_vgene,
        j_gene == target_jgene,
        cdr3_aa_length == target_cdr3_length) %>%
    dplyr::filter(!grepl("N", sequence_alignment)) %>%
    dplyr::filter(!grepl("N", sequence)) %>%
    dplyr::filter(!(startsWith(sequence_alignment,".")))

print(selected_naives)

# select on mutation frequency if provided
if (!is.null(max_mu_freq)){
    if (!("mu_freq" %in% colnames(repertoire))){
        stop("The repertoire file must contain column 'mu_freq' if you plan to select on mutation frequency.")
    }
    selected_naives <- selected_naives %>% filter(mu_freq <= max_mu_freq)
}

# select a number of sequences
if (!is.null(max_mu_freq) & "clone_id" %in% colnames(selected_naives)){
    selected_naives <- selected_naives %>%
        dplyr::group_by(clone_id, subject_id) %>%
        dplyr::slice_min(mu_freq,n=1) %>%
        dplyr::ungroup()
}

# select the number of sequences per subject
selected_naives <- selected_naives %>%
    dplyr::group_by(subject_id) %>%
    dplyr::slice_sample(n = selection_number) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(naive_sequence = gsub("\\.","",sequence_alignment)) %>%
    dplyr::mutate(seqid = paste(subject_id,sequence_id,sep="_", collapse=""))

# write to fasta file
seqs_with_id <- selected_naives$naive_sequence
names(seqs_with_id) <- selected_naives$seqid
write.fasta(sequences = as.list(seqs_with_id), names = names(seqs_with_id), file = opt$outname)
