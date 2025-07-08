process SIMULATION {
    tab "$meta.id"
    publishDir "$params.outdir/simulation/$meta.id", mode: 'copy'

    input:
    tuple val(meta), path(repertoire)

    output:
    tuple val(meta), path(simulation), path(repertoire), emit: simulation

    script:
    """
    simulation.py --input_repertoire ${repertoire} \
    --target_aa ${params.target_seq_aa} \
    --embedding_model_input ${params.embedding_model} \
    --distance_function_input ${params.distance_function} \
    --n_children_per_generation ${params.n_children} \
    --mutation_rate ${params.mutation_rate} \
    --sequence_nt_column ${params.sequence_nt_column} \
    --sequence_aa_column ${params.sequence_aa_column} \
    --output_prefix ${meta.id}
    """
}
