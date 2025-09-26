process IGBLAST {
    tag "$meta.id"
    publishDir "$params.outdir/simulation/$meta.id", mode: 'copy'
    container "community.wave.seqera.io/library/igblast:1.22.0--df7afc24896f633e"

    input:
    tuple val(meta), path(repertoire)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}_igblast.tsv"), emit: repertoire
    path "versions.yml", emit: versions

    script:
    """
    run_igblast.py -i ${repertoire} -r ${reference} -o ${meta.id}_igblast.tsv

    # Save versions
    cat <<EOF > versions.yml
    python: \$(python --version 2>&1 | head -n 1)
    pandas: \$(python -c "import pandas; print(pandas.__version__)" | head -n 1)
    EOF
    """
}
