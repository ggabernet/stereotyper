process PREPROCESS_REPERTOIRE {
    tag "$meta.id"
    label "process_low"
    publishDir "$params.outdir/simulation/$meta.id", mode: 'copy'
    container "community.wave.seqera.io/library/pandas_r-alakazam_rpy2:cea1b8fb7f0c04d3"

    input:
    tuple val(meta), path(repertoire)

    output:
    tuple val(meta), path("${meta.id}_subsampled.tsv"), emit: repertoire
    path "versions.yml", emit: versions

    script:
    """
    preprocess_repertoire.py \
    --input_repertoire ${repertoire} \
    --outname ${meta.id}_subsampled.tsv

    # Save versions
    cat <<EOF > versions.yml
    python: \$(python --version 2>&1 | head -n 1)
    pandas: \$(python -c "import pandas; print(pandas.__version__)" | head -n 1)
    r: \$(R --version | head -n 1)
    rpy2: \$(python -c "import rpy2; print(rpy2.__version__)" | head -n 1)
    EOF
    """
}
