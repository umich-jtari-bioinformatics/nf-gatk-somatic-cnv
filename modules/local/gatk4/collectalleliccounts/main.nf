nextflow.enable.dsl = 2

process GATK4_COLLECTALLELICCOUNTS {
    tag { "${meta.id}" }
    label 'COLLECT_ALLELIC_COUNTS'

    publishDir = [
        path: { "${params.outdir}/${task.process}".toString() },
        mode: params.publish_dir_mode
    ]

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(cram), path(crai), path(reference), path(reference_fai), path(reference_dict), path(sites_vcf), path(sites_vcf_tbi)

    output:
    tuple val(meta), path("${meta.id}.allelicCounts.tsv")

    when:
    cram && sites_vcf

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    gatk ${args} CollectAllelicCounts \
      -I ${cram} \
      -R ${reference} \
      -L ${sites_vcf} \
      -O ${meta.id}.allelicCounts.tsv
    """
}