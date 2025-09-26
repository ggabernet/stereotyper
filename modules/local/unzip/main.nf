process UNZIP {
    container "community.wave.seqera.io/library/unzip:6.0--0e729f0c20458893"

    input:
    path zip_file

    output:
    path "*", emit: unzipped

    script:
    """
    unzip ${zip_file}
    """
}
