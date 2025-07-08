process SELECT_SIMULATED_SEQUENCES {
    tab "$meta.id"
    publishDir "$params.outdir/simulation/$meta.id", mode: 'copy'

    input:
    tuple val(meta), path(repertoire), path(repertoire_embedding)
    tuple val(meta), path(simulated_sequences), path(simulated_embedding)


    output:
    path "${meta.id}_selected_simulated_sequences.tsv", emit: selected_simulated_sequences

    script:
    """
    select_children.py --dir . \
    --repertoire_embedding ${repertoire_embedding} \
    --simulated_embedding ${simulated_embedding} \
    --repertoire ${repertoire} \
    --simulated ${simulated_sequences} \
    --embedding_model ${params.embedding_model} \
    --target_aa ${params.target_sequence_aa} \
    --abundance ${params.abundance_fraction} \
    --repertoire_sample ${params.subsample_size} \
    --random_seed ${params.random_seed}
    """
}
