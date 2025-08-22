nextflow.enable.dsl = 2

// Interval generation (nf-core modules)
include { GATK4_PREPROCESSINTERVALS as PREPROCESS_INTERVALS } from './modules/nf-core/gatk4/preprocessintervals/main.nf'
include { GATK4_ANNOTATEINTERVALS  as ANNOTATE_INTERVALS   } from './modules/nf-core/gatk4/annotateintervals/main.nf'

// Subworkflows
include { PON_BUILD   } from './subworkflows/pon_build.nf'
include { SOMATIC_CNV } from './subworkflows/somatic_cnv.nf'

/*
 * Params (configured in nextflow.config)
 *  --samplesheet
 *  --reference_fasta
 *  --intervals (optional)
 *  --mode ('wgs'|'wes'; default in config)
 *  --bin_length (WGS auto-intervals)
 *  --padding (WES capture padding)
 *  --common_snps_vcf (optional)
 *  --build_pon_only (bool)
 *  --autosomes_only (bool)
 */
def checkParams() {
    if (!params.samplesheet)     exit 1, "ERROR: --samplesheet is required"
    if (!params.reference_fasta) exit 1, "ERROR: --reference_fasta is required"
    if (!(params.mode in ['wes','wgs'])) exit 1, "ERROR: --mode must be 'wes' or 'wgs'"
}

def parseSamplesheet(csv) {
    Channel
        .fromPath(csv)
        .splitCsv(header:true)
        .map { row ->
            // Expect: sample_id,type,cram,crai,sex,tumor_normal_id
            def sid   = row.sample_id
            def type  = row.type?.toLowerCase()
            def cram  = file(row.cram)
            def crai  = file(row.crai)
            def sex   = row.sex ?: ''
            def tnorm = row.tumor_normal_id ?: ''
            def meta  = [ id: sid, sample: sid, type: type, sex: sex ]
            tuple(meta, cram, crai, tnorm)
        }
}

workflow {
    checkParams()

    reference_fasta = file(params.reference_fasta)

    // 1) Samples
    all_samples = parseSamplesheet(params.samplesheet)

    normals_ch = all_samples
        .filter { meta, cram, crai, tnorm -> meta.type == 'normal' }
        .map { meta, cram, crai, tnorm -> tuple(meta, cram, crai) }

    tumors_ch = all_samples
        .filter { meta, cram, crai, tnorm -> meta.type == 'tumor' }
        .map { meta, tcram, tcrai, tnorm ->
            // Resolve matched normal files if provided
            def has_normal = tnorm ? true : false
            def nrec = has_normal ? all_samples
                .filter { m2, c2, i2, t2 -> m2.id == tnorm && m2.type == 'normal' }
                .first()
                .ifEmpty { null } : null

            def ncram = nrec ? nrec[1] : null
            def ncrai = nrec ? nrec[2] : null

            tuple(meta, tcram, tcrai, has_normal, ncram, ncrai)
        }

    // 2) Intervals: provided or generated (preprocess + annotate)
    Channel
        .value(params.intervals ? file(params.intervals) : null)
        .ifEmpty {
            def cap_bed = (params.mode == 'wes' && params.intervals) ? file(params.intervals) : null
            pre = PREPROCESS_INTERVALS(
                reference: reference_fasta,
                mode: params.mode,
                bin_length: params.bin_length as int,
                padding: params.padding as int,
                capture_bed: cap_bed
            )
            ann = ANNOTATE_INTERVALS(
                reference: reference_fasta,
                preprocessed_intervals: pre.out
            )
            return ann.out
        }
        .set { intervals }

    // 3) Build PoN if normals exist
    //def has_normals = normals_ch.count().tap { it }.view()
    // count normals (channel that emits a single integer)
    def normals_count_ch = normals_ch.count()

    // log the count
    normals_count_ch.view { n -> "Normal samples: ${n}" }

    // boolean-as-channel you can join/map with later
    def has_normals = normals_count_ch.map { n -> n > 0 }
    pon_ch = has_normals > 0 ? PON_BUILD(normals_ch, intervals, reference_fasta).pon_hdf5 : Channel.empty()

    if (params.build_pon_only) {
        pon_ch.view { "PON built: ${it}" }
        emit:
            pon = pon_ch
        return
    }

    // 4) Somatic CNV
    common_snps = params.common_snps_vcf ? file(params.common_snps_vcf) : null

    results = SOMATIC_CNV(tumors_ch, intervals, reference_fasta, pon_ch, common_snps)

    emit:
        segments        = results.segments
        allelicSegments = results.allelic_segments
        calls           = results.calls
}