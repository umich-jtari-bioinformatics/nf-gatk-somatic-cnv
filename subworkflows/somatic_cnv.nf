// Include GATK4 modules for the somatic CNV workflow
include { GATK4_COLLECTREADCOUNTS as GATK4_COLLECTREADCOUNTS_TUMOR }      from '../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_DENOISEREADCOUNTS }      from '../modules/nf-core/gatk4/denoisereadcounts/main.nf'
include { GATK4_COLLECTALLELICCOUNTS as GATK4_COLLECTALLELICCOUNTS_TUMOR }   from '../modules/local/gatk4/collectalleliccounts/main.nf'
include { GATK4_COLLECTALLELICCOUNTS as GATK4_COLLECTALLELICCOUNTS_NORMAL }  from '../modules/local/gatk4/collectalleliccounts/main.nf'
include { GATK4_MODELSEGMENTS }          from '../modules/nf-core/gatk4/modelsegments/main.nf'
include { GATK4_CALLCOPYRATIOSEGMENTS }  from '../modules/local/gatk4/callcopyratiosegments/main.nf'

workflow SOMATIC_CNV {

    take:
    tumors_ch         // [ meta, cram, crai, tumor_normal_id ]
    normals_ch        // [ meta, cram, crai ] - all normal samples
    pon_ch            // PoN HDF5 file channel
    intervals         // path (BED or interval_list)
    snp_vcf           // path to .vcf.gz (and .tbi alongside)
    reference_fasta   // path to reference fasta
    reference_fai     // path to reference fai
    reference_dict    // path to reference dict

    main:
    // Prepare normals for joining by sample ID
    normals_by_id = normals_ch
        .map { meta, cram, crai -> tuple(meta.id, meta, cram, crai) }

    // 1. CollectReadCounts for tumor samples
    tumor_samples = tumors_ch.map { meta, cram, crai, tnorm -> tuple(meta, cram, crai) }
    count_in = tumor_samples.combine(intervals)
    
    // Reference files as separate channels (matching pon_build.nf pattern)
    fasta_ch = Channel.value([[:], reference_fasta])
    fai_ch = Channel.value([[:], reference_fai])
    dict_ch = Channel.value([[:], reference_dict])
    
    tumor_counts = GATK4_COLLECTREADCOUNTS_TUMOR(count_in, fasta_ch, fai_ch, dict_ch)

    // 2. DenoiseReadCounts using PoN
    denoised = GATK4_DENOISEREADCOUNTS(
        tumor_counts.hdf5,
        pon_ch.map { pon -> tuple([id: 'pon'], pon) }
    )

    // 3. CollectAllelicCounts for both tumor and matched normal samples
    snp_vcf_tbi = snp_vcf.toString() + '.tbi'
    
    // Collect allelic counts for tumors
    tumor_samples_with_ref = tumors_ch.map { meta, cram, crai, tnorm -> 
        tuple(meta, cram, crai, reference_fasta, reference_fai, reference_dict, snp_vcf, file(snp_vcf_tbi)) 
    }
    
    tumor_allelic_counts = GATK4_COLLECTALLELICCOUNTS_TUMOR(
        tumor_samples_with_ref
    )
    
    // Join tumors with their matched normals for allelic counts
    tumors_with_normals = tumors_ch
        .filter { meta, cram, crai, tnorm -> tnorm && tnorm != '' }
        .map { meta, cram, crai, tnorm -> tuple(tnorm, meta, cram, crai) }
        .join(normals_by_id, by: 0)
        .map { tnorm, tumor_meta, tumor_cram, tumor_crai, normal_meta, normal_cram, normal_crai ->
            // Create new meta with tumor ID prefix to track pairing
            def paired_meta = [id: "${tumor_meta.id}_normal", sample: normal_meta.sample, type: 'normal', tumor_id: tumor_meta.id]
            tuple(paired_meta, normal_cram, normal_crai, reference_fasta, reference_fai, reference_dict, snp_vcf, file(snp_vcf_tbi))
        }
    
    matched_normal_samples = tumors_with_normals
    
    normal_allelic_counts = GATK4_COLLECTALLELICCOUNTS_NORMAL(
        matched_normal_samples
    )

    // 4. ModelSegments - prepare input combining denoised counts with allelic counts
    // Handle both tumor-only and matched tumor-normal modes
    
    // First, create a map of tumor allelic counts by sample ID
    tumor_allelic_by_id = tumor_allelic_counts
        .map { meta, allelic_tsv -> 
            println "DEBUG: Tumor allelic - ID: ${meta.id}, file: ${allelic_tsv}"
            tuple(meta.id, meta, allelic_tsv) 
        }
    
    // Create a map of normal allelic counts by tumor_id
    normal_allelic_by_tumor_id = normal_allelic_counts
        .map { meta, allelic_tsv -> 
            println "DEBUG: Normal allelic - tumor_id: ${meta.tumor_id}, file: ${allelic_tsv}"
            tuple(meta.tumor_id, meta, allelic_tsv) 
        }
    
    // Join tumor allelic counts with matched normal allelic counts (remainder: true for tumor-only samples)
    paired_allelic_counts = tumor_allelic_by_id
        .join(normal_allelic_by_tumor_id, by: 0, remainder: true)
        .map { items ->
            def tumor_id = items[0]
            def tumor_meta = items[1]
            def tumor_allelic_tsv = items[2]
            def normal_meta = items.size() > 3 ? items[3] : null
            def normal_allelic_tsv = items.size() > 4 ? items[4] : []
            
            println "DEBUG: Joined - tumor_id: ${tumor_id}, normal_allelic_tsv: ${normal_allelic_tsv}"
            tuple(tumor_meta, tumor_allelic_tsv, normal_allelic_tsv)
        }
    
    // Combine with denoised counts
    model_input = denoised.denoised
        .join(paired_allelic_counts, by: 0)
        .map { meta, denoised_tsv, tumor_allelic_tsv, normal_allelic_tsv -> 
            tuple(meta, denoised_tsv, tumor_allelic_tsv, normal_allelic_tsv)
        }

    segments = GATK4_MODELSEGMENTS(
        model_input
    )

    // 5. CallCopyRatioSegments
    copyratio_calls = GATK4_CALLCOPYRATIOSEGMENTS(
        segments.cr_seg
    )

    emit:
    segments         = segments.segmented
    copyratio        = copyratio_calls
    denoised_counts  = denoised.denoised
    tumor_allelic_counts = tumor_allelic_counts
    normal_allelic_counts = normal_allelic_counts
}