process IGBLAST {
    tag "$meta.id"
    label "process_low"
    publishDir "$params.outdir/simulation/$meta.id", mode: 'copy'
    container "community.wave.seqera.io/library/igblast_pandas:1e4fe10623287fd0"

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
    igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
    EOF
    """
}
