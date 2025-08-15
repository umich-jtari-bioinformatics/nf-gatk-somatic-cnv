nextflow.enable.dsl = 2

process GATK4_CALLCOPYRATIOSEGMENTS {
    tag { "${meta.id}" }
    label 'CALL_COPYRATIO_SEGMENTS'

    publishDir = [
        path: { "${params.outdir}/${task.process}".toString() },
        mode: params.publish_dir_mode
    ]

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(segments_file)

    output:
    tuple val(meta), path("${meta.id}.calls.cns.tsv")

    when:
    segments_file

    script:
    """
    set -euo pipefail

    gatk CallCopyRatioSegments \
      --input ${segments_file} \
      --output ${meta.id}.calls.cns.tsv
    """
}