process UNZIP {
    tag "$meta.id"
    publishDir "$params.outdir/simulation/$meta.id", mode: 'copy'
    container "docker.io/library/alpine:latest"

    input:
    path zip_file

    output:
    path "*", emit: unzipped

    script:
    """
    unzip ${zip_file}
    """
}
