process SIMULATE_CLONE {
    tag "$meta.id"
    publishDir "$params.outdir/simulate_clone/$meta.id", mode: 'copy'
    container "ggabernet/bcrphylosimulation:0.4"
    label "process_single_long"

    input:
    tuple val(meta), path(naive_seq)

    output:
    tuple val(meta), path("${meta.id}_${naive_seq.baseName}.fasta"), emit: fasta
    path "versions.yml"                      , emit: versions
    path " ${meta.id}_${naive_seq.baseName}.log"                    , emit: log
    path "${meta.id}*"

    script:
    """
    TMPDIR=/tmp xvfb-run python /bcr-phylo-benchmark/bin/simulator.py \\
    --mutability /bcr-phylo-benchmark/motifs/Mutability_S5F.csv \\
    --substitution /bcr-phylo-benchmark/motifs/Substitution_S5F.csv \\
    --outbase "${meta.id}_${naive_seq.baseName}" \\
    --lambda0 0.365 \\
    --n_to_sample 500 \\
    --naive_seq_file $naive_seq \\
    --target_seq $params.target_seq_aa \\
    --carry_cap 1000 \\
    --skip_update 100 \\
    --obs_times 100 1000 2000 3000 4000 5000 6000 \\
    --stop_dist 1 \\
    --debug 1 \\
    --multifurcating_tree \\
    --n_tries 1000 > ${meta.id}_${naive_seq.baseName}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( echo \$(python --version | grep -o "[0-9\\. ]\\+") )
        pandas: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)"))
    END_VERSIONS
    """
}
