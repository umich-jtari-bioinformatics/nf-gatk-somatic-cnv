nextflow.enable.dsl = 2

include { GATK4_COLLECTREADCOUNTS as GATK4_COLLECTREADCOUNTS_PON } from '../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_CREATEREADCOUNTPANELOFNORMALS } from '../modules/nf-core/gatk4/createreadcountpanelofnormals/main.nf'

workflow PON_BUILD {
    take:
        normals_ch       // tuple(meta, cram, crai)
        intervals        // path
        reference_fasta  // path
        reference_fai    // path
        reference_dict   // path

    main:
        // Prepare input channels for GATK4_COLLECTREADCOUNTS
        count_in = normals_ch.combine(intervals)
        
        // Reference files as separate channels with dummy meta
        fasta_ch = Channel.value([[:], reference_fasta])
        fai_ch = Channel.value([[:], reference_fai])
        dict_ch = Channel.value([[:], reference_dict])

        counts = GATK4_COLLECTREADCOUNTS_PON(count_in, fasta_ch, fai_ch, dict_ch)

        // Collect all HDF5 files and create a single tuple for PoN creation
        collected_counts = counts.hdf5
            .map { meta, hdf5 -> hdf5 }
            .collect()
            .map { hdf5_list -> tuple([id: 'pon'], hdf5_list) }

        pon = GATK4_CREATEREADCOUNTPANELOFNORMALS(collected_counts)

    emit:
        pon_hdf5 = pon.pon
}