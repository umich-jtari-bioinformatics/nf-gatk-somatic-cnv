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
    Channel
      .of( params.intervals ? file(params.intervals) : null )
      .map { it ->
          if( it == null ) {
              // No intervals path given; for WGS you could generate by binning, but we enforce an explicit file
              exit 1, "ERROR: --intervals is required (BED or interval_list)."
          }
          def p = it as File
          def name = p.getName().toLowerCase()
          if( name.endsWith('.interval_list') || name.endsWith('.intervals') ) {
              // Already interval_list → pass through
              tuple(p, null) // (interval_list, annotated)
          } else if( name.endsWith('.bed') ) {
              // WES BED → preprocess + annotate via nf-core modules
              def pre = PREPROCESS_INTERVALS(
                  reference: reference_fasta,
                  intervals: p,
                  // module expects specific flags; bin_length is ignored for WES, padding used
                  mode: params.mode ?: 'wes',
                  bin_length: (params.bin_length ?: 1000) as int,
                  padding: (params.padding ?: 0) as int
              ).out

              def ann = ANNOTATE_INTERVALS(
                  reference: reference_fasta,
                  preprocessed_intervals: pre
              ).out

              tuple(pre, ann)
          } else {
              exit 1, "ERROR: --intervals must be .interval_list/.intervals or .bed, got: ${p}"
          }
      }
      .set { interval_pair }

    intervals_ch = interval_pair.map { pre, ann -> pre }

    // 3) Build PoN from normals

    pon_ch = PON_BUILD(
        normals_ch,
        intervals_ch,
        reference_fasta
    )

    if( params.build_pon_only as boolean ) {
        pon_ch.view { p -> "PON written: ${p}" }
        emit:
            pon = pon_ch
        return
    }

    // 4) Somatic CNV on tumors
    def common_snps = params.common_snps_vcf ? file(params.common_snps_vcf) : null

    somatic = SOMATIC_CNV(
        tumors_ch: tumors,          // (sid, cram, crai, tnid)
        normals_map: normals_map,   // (nid, [cram, crai])  <-- channel
        fasta: fasta, dict: dict, fai: fai,
        intervals: intervals_ch,
        pon_hdf5: pon,
        snp_vcf: file(params.snp_vcf)  // or null if skipping
    )

    emit:
        cr_segments      = somatic.cr_segments
        calls            = somatic.calls
}
