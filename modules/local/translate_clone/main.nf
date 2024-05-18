process TRANSLATE_CLONE {
    tag "${meta.id}_${meta.naive}"
    publishDir "$params.outdir/add_clone_repertoire/$meta.id", mode: 'copy'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:a9ee25632c9b10bbb012da76e6eb539acca8f9cd-1' :
        'quay.io/biocontainers/mulled-v2-7d8e418eb73acc6a80daea8e111c94cf19a4ecfd:a9ee25632c9b10bbb012da76e6eb539acca8f9cd-1' }"
    label "process_medium"

    input:
    tuple val(meta), path(clone)

    output:
    tuple val(meta), path("${meta.id}_${meta.naive}_translated.tsv") , emit: airr

    script:
    """
    wget -c https://github.com/nf-core/test-datasets/raw/airrflow/database-cache/igblast_base.zip
    unzip igblast_base.zip

    igblastn -germline_db_V igblast_base/database/imgt_human_ig_v \\
    -germline_db_D igblast_base/database/imgt_human_ig_d \\
    -germline_db_J igblast_base/database/imgt_human_ig_j \\
    -auxiliary_data igblast_base/optional_file/human_gl.aux \\
    -organism human \\
    -show_translation \\
    -query ${clone} \\
    -outfmt 19 \\
    -out ${meta.id}_${meta.naive}_translated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igblastn: \$( igblastn -version | grep -o "igblast[0-9\\. ]\\+" | grep -o "[0-9\\. ]\\+" )
    END_VERSIONS
    """


}
