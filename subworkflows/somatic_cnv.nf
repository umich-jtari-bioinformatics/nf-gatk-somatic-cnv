workflow SOMATIC_CNV {

    take:
    tumors_ch         // [ meta, cram, crai ]
    normals_ch_opt    // [ meta, cram, crai ] or Channel.empty()
    intervals         // path (BED or interval_list)
    snp_vcf           // path to .vcf.gz (and .tbi alongside)

    main:
    // placeholder channels so the file parses before you wire modules
    segments_ch      = Channel.empty()
    copyratio_calls  = Channel.empty()

    emit:
    segments         = segments_ch
    copyratio        = copyratio_calls
}