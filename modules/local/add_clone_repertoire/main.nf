process ADD_CLONE_REPERTOIRE {
    tag "${meta.id}_${meta.naive}_${quant}"
    publishDir "$params.outdir/repertoire_with_clones/${meta.id}/"
    container "immcantation/airrflow:3.2.0"
    label "process_single"

    input:
    tuple val(quant), val(meta), path(repertoire), path(clone)

    output:
    tuple val(meta), path("${meta.id}_repertoire_with_${meta.naive}_clone.tsv")

    script:
    """
    addclone_to_rep.R --repertoire ${repertoire} \\
    --clone ${clone} \\
    --quant ${quant} \\
    --abund $params.clonal_abundance \\
    --target_sequence $params.target_seq_aa \\
    --outname ${meta.id}_repertoire_with_${meta.naive}_clone_sampled_${quant}.tsv
    """
}
