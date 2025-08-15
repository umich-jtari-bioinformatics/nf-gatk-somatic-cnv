nextflow.enable.dsl = 2

include { GATK4_COLLECTREADCOUNTS     } from '../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_DENOISEREADCOUNTS     } from '../modules/nf-core/gatk4/denoisereadcounts/main.nf'
include { GATK4_MODELSEGMENTS         } from '../modules/nf-core/gatk4/modelsegments/main.nf'

// Local modules for the two missing tools
include { GATK4_COLLECTALLELICCOUNTS  } from '../modules/local/gatk4/collectalleliccounts/main.nf'
include { GATK4_CALLCOPYRATIOSEGMENTS } from '../modules/local/gatk4/callcopyratiosegments/main.nf'

workflow SOMATIC_CNV {
    take:
        tumors_ch        // tuple(meta, tcram, tcrai, has_normal, ncram, ncrai)
        intervals        // path
        reference_fasta  // path
        pon_hdf5         // path
        common_snps_vcf  // path or null

    main:
        // 1) Tumor read counts
        t_counts = GATK4_COLLECTREADCOUNTS(
            tumors_ch.map { meta, tcram, tcrai, has_normal, ncram, ncrai ->
                tuple(meta, tcram, tcrai, intervals, reference_fasta)
            }
        )

        // 2) Denoise with PoN
        t_denoised = GATK4_DENOISEREADCOUNTS(
            t_counts.out.map { meta, counts_hdf5 ->
                tuple(meta, counts_hdf5, pon_hdf5)
            }
        )

        // 3) Allelic counts (optional if common_snps_vcf given)
        t_allelic = common_snps_vcf ?
            GATK4_COLLECTALLELICCOUNTS(
                tumors_ch.map { meta, tcram, tcrai, has_normal, ncram, ncrai ->
                    tuple(meta, tcram, tcrai, reference_fasta, common_snps_vcf)
                }
            )
            : Channel.empty()

        n_allelic = (common_snps_vcf ?
            GATK4_COLLECTALLELICCOUNTS(
                tumors_ch
                    .filter { meta, tcram, tcrai, has_normal, ncram, ncrai -> has_normal && ncram != null && ncrai != null }
                    .map    { meta, tcram, tcrai, has_normal, ncram, ncrai ->
                        // keep keying by tumor meta.id; we don't change id
                        tuple(meta, ncram, ncrai, reference_fasta, common_snps_vcf)
                    }
            )
            : Channel.empty())

        // 4) Join denoised CR with allelic counts by tumor id
        den_k = t_denoised.out.map { meta, denoisedCR, standardizedCR ->
            tuple(meta.id, meta, denoisedCR)
        }
        t_ac_k = t_allelic ? t_allelic.out.map { meta, tsv -> tuple(meta.id, tsv) } : Channel.empty()
        n_ac_k = n_allelic ? n_allelic.out.map { meta, tsv -> tuple(meta.id, tsv) } : Channel.empty()

        merged = den_k
            .combine(t_ac_k, remainder: true)
            .combine(n_ac_k, remainder: true)

        // 5) ModelSegments
        modeled = GATK4_MODELSEGMENTS(
            merged.map { id, meta_den, tAllelic, nAllelic ->
                def meta  = meta_den[1]
                def denCR = meta_den[2]
                // common nf-core signature: tuple(meta, denoisedCR, tumor_allelic?, normal_allelic?, intervals)
                tuple(meta, denCR, tAllelic, nAllelic, intervals)
            }
        )

        // 6) CallCopyRatioSegments (local)
        called = GATK4_CALLCOPYRATIOSEGMENTS(
            modeled.out.map { meta, segments_file, allelic_seg_file, model_dir ->
                tuple(meta, segments_file)
            }
        )

    emit:
        segments         = modeled.out.map { meta, segments_file, allelic_seg_file, model_dir -> tuple(meta, segments_file) }
        allelic_segments = modeled.out.map { meta, segments_file, allelic_seg_file, model_dir -> tuple(meta, allelic_seg_file) }
        calls            = called.out.map  { meta, calls_file -> tuple(meta, calls_file) }
}