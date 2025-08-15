nextflow.enable.dsl = 2

include { GATK4_COLLECTREADCOUNTS             } from '../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_CREATEREADCOUNTPANELOFNORMALS } from '../modules/nf-core/gatk4/createreadcountpanelofnormals/main.nf'

workflow PON_BUILD {
    take:
        normals_ch       // tuple(meta, cram, crai)
        intervals        // path
        reference_fasta  // path

    main:
        count_in = normals_ch.map { meta, cram, crai ->
            tuple(meta, cram, crai, intervals, reference_fasta)
        }

        counts = GATK4_COLLECTREADCOUNTS(count_in)

        pon = GATK4_CREATEREADCOUNTPANELOFNORMALS(
            counts.out.map { meta, hdf5 -> hdf5 }.toList()
        )

    emit:
        pon_hdf5 = pon.out
}