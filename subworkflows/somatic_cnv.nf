/*
 * SOMATIC_CNV subworkflow
 *
 * Inputs:
 *   tumors_ch    : Channel of tuples (sid, tumor_cram, tumor_crai, tumor_normal_id)
 *   normals_map  : Channel of tuples (normal_id, [normal_cram, normal_crai])
 *   fasta, dict, fai : reference files
 *   intervals    : interval_list channel
 *   pon_hdf5     : PoN path/channel
 *   snp_vcf      : VCF path (or null when params.skip_allelic_counts is true)
 *
 * Emits:
 *   cr_segments  : tuples (sid, <path to .cr.seg>)
 *   calls        : <path to called segments>
 */

// nf-core modules: real + stub aliases
include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS } from '../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS_STUB } from '../modules/nf-core/gatk4/collectreadcounts/main.nf' stub: true

include { GATK4_DENOISEREADCOUNTS as DENOISEREADCOUNTS } from '../modules/nf-core/gatk4/denoisereadcounts/main.nf'
include { GATK4_DENOISEREADCOUNTS as DENOISEREADCOUNTS_STUB } from '../modules/nf-core/gatk4/denoisereadcounts/main.nf' stub: true

include { GATK4_COLLECTALLELICCOUNTS as COLLECTALLELICCOUNTS } from '../modules/nf-core/gatk4/collectalleliccounts/main.nf'
include { GATK4_COLLECTALLELICCOUNTS as COLLECTALLELICCOUNTS_STUB } from '../modules/nf-core/gatk4/collectalleliccounts/main.nf' stub: true

include { GATK4_MODELSEGMENTS as MODELSEGMENTS } from '../modules/nf-core/gatk4/modelsegments/main.nf'
include { GATK4_MODELSEGMENTS as MODELSEGMENTS_STUB } from '../modules/nf-core/gatk4/modelsegments/main.nf' stub: true

include { GATK4_CALLCOPYRATIOSEGMENTS as CALLCOPYRATIOSEGMENTS } from '../modules/nf-core/gatk4/callcopyratiosegments/main.nf'
include { GATK4_CALLCOPYRATIOSEGMENTS as CALLCOPYRATIOSEGMENTS_STUB } from '../modules/nf-core/gatk4/callcopyratiosegments/main.nf' stub: true


workflow SOMATIC_CNV {

    take:
    tumors_ch
    normals_map             // (normal_id, [cram, crai])
    fasta
    dict
    fai
    intervals
    pon_hdf5
    snp_vcf

    main:
    /*
     * 1) Tumor read counts
     */
    tumor_counts = tumors_ch.map { sid, tcram, tcrai, tnid ->
        def proc = params.stub ? COLLECTREADCOUNTS_STUB : COLLECTREADCOUNTS
        proc(
            sample_id: sid,
            cram: file(tcram),
            crai: file(tcrai),
            fasta: fasta, dict: dict, fai: fai,
            intervals: intervals
        ).counts.map { tuple(sid, it, tnid) }
    }.flatten()

    /*
     * 2) Denoise with PoN
     */
    denoised = tumor_counts.map { sid, counts_hdf5, tnid ->
        def proc = params.stub ? DENOISEREADCOUNTS_STUB : DENOISEREADCOUNTS
        proc(
            sample_id: sid,
            counts_hdf5: counts_hdf5,
            pon_hdf5: pon_hdf5
        ).denoised.map { tuple(sid, it, tnid) }
    }.flatten()

    /*
     * 3) Allelic counts (tumor ± matched normal)
     */
    // Tumor AC (or null when skipping)
    tumor_ac = tumors_ch.map { sid, tcram, tcrai, tnid ->
        if (params.skip_allelic_counts) {
            return tuple(sid, null)
        }
        def proc = params.stub ? COLLECTALLELICCOUNTS_STUB : COLLECTALLELICCOUNTS
        proc(
            sample_id: sid,
            cram: file(tcram),
            crai: file(tcrai),
            fasta: fasta, dict: dict, fai: fai,
            snp_vcf: snp_vcf
        ).ac.map { tuple(sid, it) }
    }.flatten()

    // Prepare a lookup channel for normals: (normal_id, ncram, ncrai)
    normals_lookup = normals_map.map { nid, pair -> tuple(nid, file(pair[0]), file(pair[1])) }

    // Key tumors by their tumor_normal_id for join: (tnid, sid, tcram, tcrai)
    tumors_keyed = tumors_ch
        .filter { sid, tcram, tcrai, tnid -> tnid != null && (tnid as String) }
        .map    { sid, tcram, tcrai, tnid -> tuple(tnid, sid, file(tcram), file(tcrai)) }

    // Inner-join tumors that have tnid with normals_lookup: (tnid, sid, tcram, tcrai, ncram, ncrai)
    tumor_normal_pairs = tumors_keyed.join(normals_lookup)

    // Normal AC (only for tumors with a tnid & when not skipping)
    normal_ac = Channel.empty()
    if( !params.skip_allelic_counts ) {
        normal_ac = tumor_normal_pairs.map { tnid, sid, tcram, tcrai, ncram, ncrai ->
            def proc = params.stub ? COLLECTALLELICCOUNTS_STUB : COLLECTALLELICCOUNTS
            proc(
                sample_id: tnid,     // name AC by the normal's ID
                cram: ncram,
                crai: ncrai,
                fasta: fasta, dict: dict, fai: fai,
                snp_vcf: snp_vcf
            ).ac.map { tuple(sid, it) }   // re-key on tumor sid for later join
        }.flatten()
    }

    /*
     * 4) Model segments
     */
    modeled = denoised
        .join(tumor_ac)        // by sid
        .map { sid1, denoised_cr, sid2, tumor_ac_file -> tuple(sid1, denoised_cr, tumor_ac_file) }
        .leftJoin(normal_ac)   // optional
        .map { sid, denoised_cr, tumor_ac_file, _sidn, normal_ac_file ->
            def proc = params.stub ? MODELSEGMENTS_STUB : MODELSEGMENTS
            proc(
                sample_id: sid,
                denoised_cr: denoised_cr,
                tumor_ac: tumor_ac_file,
                normal_ac_optional: normal_ac_file
            ).cr_seg.map { tuple(sid, it) }
        }.flatten()

    /*
     * 5) Call copy ratio segments
     */
    calls = modeled.map { sid, cr_seg ->
        def proc = params.stub ? CALLCOPYRATIOSEGMENTS_STUB : CALLCOPYRATIOSEGMENTS
        proc(sample_id: sid, cr_seg: cr_seg).calls
    }.flatten()

    emit:
    cr_segments = modeled
    calls       = calls
}