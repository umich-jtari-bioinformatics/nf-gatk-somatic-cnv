workflow SOMATIC_CNV {
    take:
    tumors_ch
    normals_ch_opt
    intervals
    snp_vcf
    ref_fasta
    ref_fai
    ref_dict
    pon_hdf5

    main:

    emit:
}
