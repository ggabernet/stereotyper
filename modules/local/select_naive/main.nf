process SELECT_NAIVES {
    tag "$meta.id"
    publishDir "$params.outdir/select_naives/$meta.id", mode: 'copy'
    container "immcantation/airrflow:3.2.0"
    label "process_single"

    input:
    tuple val(meta), path(repertoire)

    output:
    tuple val(meta), path("naive_seqs.fasta"), emit: fasta
    path "versions.yml"                      , emit: versions


    script:
    """
    select_naives.R --repertoire $repertoire \\
                    --target $params.target_seq_aa \\
                    --target_vgene $params.target_vgene \\
                    --target_jgene $params.target_jgene \\
                    --target_cdr3 $params.target_cdr3 \\
                    --max_mu_freq $params.max_mu_freq \\
                    --selection_number $params.selection_number \\
                    --outname naive_seqs.fasta

    echo "${task.process}": > versions.yml
    Rscript -e "cat(paste0('  alakazam: ',packageVersion('alakazam'),'\n'))" >> versions.yml
    Rscript -e "cat(paste0('  dplyr: ',packageVersion('dplyr'),'\n'))" >> versions.yml
    """
}
