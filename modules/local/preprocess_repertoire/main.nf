process PREPROCESS_REPERTOIRE {
    tab "$meta.id"
    publishDir "$params.outdir/simulation/$meta.id", mode: 'copy'

    input:
    tuple val(meta), path(repertoire)

    output:
    path "${meta.id}_subsampled.tsv", emit: subsampled_repertoire

    script:
    """
    python preprocess_repertoire.py \
    --input_repertoire ${repertoire} \
    --subsample_size ${params.subsample_size} \
    --outname ${meta.id}_subsampled.tsv
    """
}
