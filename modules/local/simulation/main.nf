process SIMULATE_CONVERGENCE {
    tag "$meta.id"
    publishDir "$params.outdir/simulation/$meta.id", mode: 'copy'

    input:
    tuple val(meta), path(repertoire)

    output:
    tuple val(meta), path(simulation), emit: sequences
    path "versions.yml", emit: versions

    script:
    """
    simulation.py --input_repertoire ${repertoire} \
    --target_aa ${params.target_seq_aa} \
    --distance_function_input ${params.distance_function} \
    --n_children_per_generation ${params.n_children} \
    --mutation_rate ${params.mutation_rate} \
    --sequence_nt_column ${params.sequence_nt_column} \
    --sequence_aa_column ${params.sequence_aa_column} \
    --output_prefix ${meta.id}

    # Save versions
    cat <<EOF > versions.yml
    python: \$(python --version 2>&1 | head -n 1)
    pandas: \$(python -c "import pandas; print(pandas.__version__)" | head -n 1)
    r: \$(R --version | head -n 1)
    rpy2: \$(python -c "import rpy2; print(rpy2.__version__)" | head -n 1)
    EOF
    """
}
