process MULTIQC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1bef8af6be88c5733461959c46ac8ef73d18f65277f62a1695d0e1633054f9c2/data'
        : 'community.wave.seqera.io/library/multiqc:1.34--db7c73dae76bc9e6'}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.21--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(multiqc_files, stageAs: "?/*"), path(multiqc_config, stageAs: "?/*"), path(multiqc_logo), path(replace_names), path(sample_names)

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*_data"), emit: data
    tuple val(meta), path("*_plots"), emit: plots, optional: true
    // MultiQC should not push its versions to the `versions` topic. Its input depends on the versions topic to be resolved thus outputting to the topic will let the pipeline hang forever
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${prefix} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_plots/.stub
    touch multiqc_report.html
    """
}
