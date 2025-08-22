// subworkflows/somatic_cnv.nf  (no dsl directive here)

// --- module includes (directory path, not /main.nf) ---
include { GATK4_COLLECTREADCOUNTS }     from '../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DENOISEREADCOUNTS }     from '../modules/nf-core/gatk4/denoisereadcounts/main'
include { GATK4_COLLECTALLELICCOUNTS }  from '../modules/local/gatk4/collectalleliccounts/main'
include { GATK4_MODELSEGMENTS }         from '../modules/nf-core/gatk4/modelsegments/main'
include { GATK4_CALLCOPYRATIOSEGMENTS } from '../modules/local/gatk4/callcopyratiosegments/main'

workflow SOMATIC_CNV {

  take:
    tumors_ch        // [meta, cram, crai]
    normals_ch_opt   // [meta, cram, crai] or Channel.empty()
    intervals        // BED/interval_list
    snp_vcf          // .vcf.gz (+ .tbi)
    ref_fasta
    ref_fai
    ref_dict
    pon_hdf5         // HDF5 PoN (Channel.empty() if not used)

  main:

    // Tumor read counts -> denoise with PoN
    t_counts = GATK4_COLLECTREADCOUNTS(
      tumors_ch.map { meta, cram, crai ->
        tuple(meta, cram, crai, intervals, ref_fasta, ref_fai, ref_dict)
      }
    )

    t_denoised = GATK4_DENOISEREADCOUNTS(
      t_counts.out.hdf5.map { meta, hdf5 -> tuple(meta, hdf5, pon_hdf5) }
    )
    // Adjust field names if your module uses "standardized" vs "standardised"
    t_scr_ch = t_denoised.out.standardised
    t_dcr_ch = t_denoised.out.denoised
    t_rc_h5  = t_counts.out.hdf5

    // Tumor allelic counts
    t_allelic = GATK4_COLLECTALLELICCOUNTS(
      tumors_ch.map { meta, cram, crai ->
        tuple(meta, cram, crai, snp_vcf, ref_fasta, ref_fai, ref_dict)
      }
    )
    t_ac_ch = t_allelic.out.allelic

    // Optional normal allelic counts (kept as a separate keyed channel)
    n_ac_in  = normals_ch_opt.map { meta, cram, crai ->
                  tuple(meta, cram, crai, snp_vcf, ref_fasta, ref_fai, ref_dict)
               }
    n_allelic = GATK4_COLLECTALLELICCOUNTS(n_ac_in)
    n_ac_kv   = n_allelic.out.allelic.map { meta, ac -> tuple(meta.id, ac) }   // (key=id, val=ac)

    // Keyed join of tumor SCR + tumor AC by sample id
    t_scr_kv = t_scr_ch.map { meta, scr -> tuple(meta.id, tuple(meta, scr)) }  // (id, (meta,scr))
    t_ac_kv  = t_ac_ch .map { meta, ac  -> tuple(meta.id, ac) }                // (id, ac)

    joined   = t_scr_kv.join(t_ac_kv).map { id, meta_scr, ac ->
                 def (meta, scr) = meta_scr
                 tuple(meta, scr, ac)                                          // (meta, scr, tac)
               }

    // Two parallel model paths: with and without normal AC
    model_in_no_n  = joined.map { meta, scr, tac -> tuple(meta, scr, tac) }
    modeled_no_n   = GATK4_MODELSEGMENTS(model_in_no_n)

    // If normals present, join by id; if empty, this produces no items (that’s fine)
    model_in_with_n = joined.join(n_ac_kv).map { id, (meta, scr, tac), nac ->
                        tuple(meta, scr, tac, nac)
                      }
    modeled_with_n  = GATK4_MODELSEGMENTS(model_in_with_n)

    // Merge whichever path produced results (normals path yields items only if normals exist)
    segments_ch = modeled_with_n.out.segments.mix(modeled_no_n.out.segments)

    // Call copy-ratio segments
    copyratio_calls = GATK4_CALLCOPYRATIOSEGMENTS(
      modeled_no_n.out.copy_ratio_segments.mix(modeled_with_n.out.copy_ratio_segments)
    ).out.copy_ratio_calls

  emit:
    segments   = segments_ch
    copyratio  = copyratio_calls
    readcounts = t_rc_h5
    standardised = t_scr_ch
    denoised   = t_dcr_ch
    tumor_allelic = t_ac_ch
}