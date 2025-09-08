nextflow.enable.dsl=2

// Interval generation (nf-core modules)
include { GATK4_PREPROCESSINTERVALS as PREPROCESS_INTERVALS } from "${projectDir}/modules/nf-core/gatk4/preprocessintervals/main.nf"
include { GATK4_ANNOTATEINTERVALS  as ANNOTATE_INTERVALS   } from "${projectDir}/modules/nf-core/gatk4/annotateintervals/main.nf"

// GATK4 modules used by subworkflows (import here to avoid parser grief)
include { GATK4_COLLECTREADCOUNTS }      from "${projectDir}/modules/nf-core/gatk4/collectreadcounts/main.nf"
include { GATK4_DENOISEREADCOUNTS }      from "${projectDir}/modules/nf-core/gatk4/denoisereadcounts/main.nf"
include { GATK4_COLLECTALLELICCOUNTS }   from "${projectDir}/modules/local/gatk4/collectalleliccounts/main.nf"
include { GATK4_MODELSEGMENTS }          from "${projectDir}/modules/nf-core/gatk4/modelsegments/main.nf"
include { GATK4_CALLCOPYRATIOSEGMENTS }  from "${projectDir}/modules/local/gatk4/callcopyratiosegments/main.nf"

// Subworkflows from this repo
include { PON_BUILD }  from "${projectDir}/subworkflows/pon_build.nf"
include { SOMATIC_CNV } from "${projectDir}/subworkflows/somatic_cnv.nf"

/*
 * Expected params (set in nextflow.config):
 *   --samplesheet                CSV with sample_id,type,cram,crai,sex,tumor_normal_id
 *   --reference_fasta            path to .fa/.fasta (with .fai/.dict provided separately in config to modules)
 *   --reference_dict             path to .dict
 *   --reference_fai              path to .fai
 *   --intervals                  EITHER a .interval_list/.intervals OR a .bed (BED triggers preprocess+annotate)
 *   --mode                       'wes' or 'wgs' (used only when intervals are BED or absent)
 *   --bin_length                 WGS bin size (when generating)
 *   --padding                    WES padding (when generating from BED)
 *   --common_snps_vcf            VCF for allelic counts (optional if skipping allelic counts in your config)
 *   --build_pon_only             true to stop after PON creation (matches your config)
 *   --pon_hdf5                   Pre-built panel of normals HDF5 file (optional)
 *   --outdir                     output directory
 */

def checkParams() {
    if( !params.samplesheet )         exit 1, "ERROR: --samplesheet is required."
    if( !params.reference_fasta )     exit 1, "ERROR: --reference_fasta is required."
    if( !params.reference_fai )       exit 1, "ERROR: --reference_fai is required."
    if( !params.reference_dict )      exit 1, "ERROR: --reference_dict is required."
    if( !params.intervals && params.mode == 'wes' ) {
        exit 1, "ERROR: WES mode requires --intervals pointing to a BED or an interval_list."
    }
    if( !params.pon_hdf5 && !params.build_pon_only ) {
        log.warn "No pre-built PoN provided. Will need normal samples to build PoN."
    }
}

def parseSamplesheet(csv) {
    Channel
        .fromPath(csv)
        .splitCsv(header:true)
        .map { row ->
            // Expect: sample_id,type,cram,crai,sex,tumor_normal_id
            def sid   = row.sample_id?.toString()
            def type  = row.type?.toString()?.toLowerCase()
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

    // Resolve outdir to absolute path relative to launch directory
    params.outdir = file(params.outdir).toAbsolutePath()

    // Refs
    def reference_fasta = file(params.reference_fasta)
    def reference_dict  = file(params.reference_dict)
    def reference_fai   = file(params.reference_fai)

    // 1) Samples → channels
    all_samples = parseSamplesheet(params.samplesheet)

    normals_ch = all_samples
        .filter { meta, cram, crai, tnorm -> meta.type == 'normal' }
        .map    { meta, cram, crai, tnorm -> tuple(meta, cram, crai) }

    tumors_ch = all_samples
        .filter { meta, cram, crai, tnorm -> meta.type == 'tumor' }
        .map    { meta, cram, crai, tnorm -> tuple(meta, cram, crai, tnorm) }

    // Map of normals by sample_id to find matched normals later
    normals_map = normals_ch
        .map { meta, cram, crai -> tuple(meta instanceof Map ? meta.id : meta, [cram, crai]) }

    // 2) Intervals: use provided interval_list directly, or preprocess+annotate if BED
    def intervals_file = params.intervals ? file(params.intervals) : null
    if( !intervals_file ) {
        exit 1, "ERROR: --intervals is required (BED or interval_list)."
    }
    
    def name = intervals_file.getName().toLowerCase()
    if( name.endsWith('.interval_list') || name.endsWith('.intervals') ) {
        // Already interval_list → use directly
        intervals_ch = Channel.value(intervals_file)
    } else if( name.endsWith('.bed') ) {
        // BED → preprocess + annotate via nf-core modules
        fasta_ch = Channel.value(tuple([id: 'ref'], reference_fasta))
        fai_ch = Channel.value(tuple([[:], reference_fai]))
        dict_ch = Channel.value(tuple([[:], reference_dict]))
        bed_ch = Channel.value(tuple([[:], intervals_file]))
        exclude_ch = Channel.value(tuple([[:], []]))  // empty exclude intervals
        
        preprocessed = PREPROCESS_INTERVALS(
            fasta_ch,
            fai_ch,
            dict_ch,
            bed_ch,
            exclude_ch
        )
        
        intervals_ch = preprocessed.interval_list.map { meta, interval_list -> interval_list }
    } else {
        exit 1, "ERROR: --intervals must be .interval_list/.intervals or .bed, got: ${intervals_file}"
    }

    // 3) PoN handling: use pre-built or build new
    if( params.pon_hdf5 ) {
        // Use pre-built PoN
        pon_ch = Channel.value(file(params.pon_hdf5))
        log.info "Using pre-built PoN: ${params.pon_hdf5}"
    } else {
        // Build PoN from normals
        pon_ch = PON_BUILD(
            normals_ch,
            intervals_ch,
            reference_fasta,
            reference_fai,
            reference_dict
        )
    }

    if( params.build_pon_only as boolean ) {
        pon_ch.view { p -> "PON written: ${p}" }
    } else {
        // 4) Somatic CNV on tumors
        def common_snps = params.common_snps_vcf ? file(params.common_snps_vcf) : null

        somatic = SOMATIC_CNV(
            tumors_ch,
            normals_ch,       // Pass normals for matched analysis
            pon_ch,           // Pass PoN to somatic workflow
            intervals_ch,
            file(params.snp_vcf),  // snp_vcf
            reference_fasta,
            reference_fai,
            reference_dict
        )
    }
}
