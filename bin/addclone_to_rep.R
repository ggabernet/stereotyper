#!/usr/bin/env Rscript
# Written by Gisela Gabernet and released under the MIT license (2020).


# Select naive sequences from a repertoire:
# Arguments:
#   --repertoire    Tabulated data in AIRR (TSV) format with clonal assignments and germline assignments.
#   --clone         Translated clone sequences to add to repertoire
#   --quant           Hamming distance percentile
#   --target_sequence      Amino acid sequence to be used as target for the convergence process.
#   --abund         Clonal abundance of added clone in repertoire.
#   --outname       Output file name.
#   -h  Display help.
# Example: ./addclone_to_rep.R --repertoire AIRR_repertoire.tsv --clone airr_translated_clone.tsv --quant q1 --abund 1 --outname "repertoire_with_clone.tsv"

# Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(airr))

# Define commmandline arguments
opt_list <- list(
    make_option(c("--repertoire"), default=NULL,
                help="Input repertoire .tsv file after clone definition and germline definition."),
    make_option(c("--clone"), default=NULL,
                help="Translated clone sequences to add to repertoire."),
    make_option(c("--target_sequence"), default=NULL,
                help="Amino acid sequence to be used as target for the convergence process."),
    make_option(c("--abund"), default=NULL,
                help="Clonal abundance of added clone in repertoire."),
    make_option(c("--outname"), help="Output file name.")
)
opt <- parse_args(OptionParser(option_list=opt_list))

# Read repertoire
repertoire <- airr::read_rearrangement(opt$repertoire)

# Read target sequence
target_sequence <- opt$target_sequence

# Read clone
clone <- airr::read_rearrangement(opt$clone) %>%
            dplyr::filter(productive == TRUE) %>%
            dplyr::filter(vj_in_frame == TRUE) %>%
            dplyr::filter(nchar(sequence_aa) == nchar(target_sequence))

# abundance
abund <- as.numeric(opt$abund)
n_sample <- round(nrow(repertoire) * abund, 0)

print(paste("Number of sequences in repertoire:", nrow(repertoire)))

# Getting distance of clone sequence_aa to target sequence
clone$dist_to_target <- sapply(clone$sequence_aa, alakazam::seqDist, seq2=target_sequence,
                                        dist_mat = getAAMatrix(gap=0))

# Getting quantiles
quantiles <- quantile(clone$dist_to_target, probs = c(0.25, 0.5, 0.75))
print(quantiles)

n_sample <- min(n_sample, nrow(clone))

for (i in length(quantiles)) {
        clone$sampling <- ifelse(clone$dist_to_target <= quantiles[i], 1, 0)
        sampled_clone <- clone %>% filter(sampling == 1) %>% sample_n(n_sample)
        print(paste("Number of sequences sampled for quantile ",i,":", nrow(sampled_clone)))
        repertoire_with_clone <- dplyr::bind_rows(repertoire, sampled_clone)
        write.table(repertoire_with_clone, file = paste0(opt$outname,"_quant",i,".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Sample closest
sampled_clone <- clone %>% slice_min(dist_to_target, n = n_sample)
print(paste("Number of sequences sampled for closest:", nrow(sampled_clone)))
repertoire_with_clone <- dplyr::bind_rows(repertoire, sampled_clone)
write.table(repertoire_with_clone, file = paste0(opt$outname,"_closest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Sample random
sampled_clone <- clone %>% sample_n(n_sample)
print(paste("Number of sequences sampled for random:", nrow(sampled_clone)))
repertoire_with_clone <- dplyr::bind_rows(repertoire, sampled_clone)
write.table(repertoire_with_clone, file = paste0(opt$outname,"_random.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

