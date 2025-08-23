# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Nextflow DSL2 workflow implementing the GATK somatic CNV (Copy Number Variation) pipeline. The workflow follows this sequence:
`CollectReadCounts → CreateReadCountPanelOfNormals (PoN) → DenoiseReadCounts → CollectAllelicCounts → ModelSegments → CallCopyRatioSegments`

## Common Commands

### Running the Workflow

**Build Panel of Normals (PoN) only:**
```bash
nextflow run main.nf \
  -profile docker \
  --samplesheet samples/samplesheet_normals.csv \
  --reference_fasta /path/to/reference.fasta \
  --reference_fai /path/to/reference.fasta.fai \
  --reference_dict /path/to/reference.dict \
  --intervals /path/to/intervals.bed \
  --snp_vcf /path/to/gnomad.snps.vcf.gz \
  --outdir ./results_pon \
  --build_pon_only true
```

**Full somatic CNV analysis (PoN + tumor analysis):**
```bash
nextflow run main.nf \
  -profile docker \
  --samplesheet samples/samplesheet_cohort.csv \
  --reference_fasta /path/to/reference.fasta \
  --reference_fai /path/to/reference.fasta.fai \
  --reference_dict /path/to/reference.dict \
  --intervals /path/to/intervals.bed \
  --snp_vcf /path/to/gnomad.snps.vcf.gz \
  --outdir ./results
```

**Assay-specific flags:**
- WES without intervals: add `--assay wes --capture_bed /path/capture.bed` (omit `--intervals`)
- WGS without intervals: add `--assay wgs --bin_length 1000 --interval_padding 0` (omit `--intervals`)

### Available Profiles
- `standard`: Local execution (default)
- `docker`: Docker container execution (broadinstitute/gatk:4.5.0.0)
- `singularity`: Singularity container execution
- `conda`: Conda environment (gatk4=4.5.0.0)
- `test`: Test configuration

## Architecture Overview

### Main Workflow (`main.nf`)
The main workflow orchestrates two key subworkflows:
1. **PON_BUILD**: Creates Panel of Normals from normal samples
2. **SOMATIC_CNV**: Processes tumor samples using the PoN (currently a stub/placeholder)

### Key Components

**Modules:**
- `modules/nf-core/gatk4/`: Standard nf-core GATK4 modules (collectreadcounts, createreadcountpanelofnormals, denoisereadcounts, modelsegments, preprocessintervals, annotateintervals)
- `modules/local/gatk4/`: Custom local modules (collectalleliccounts, callcopyratiosegments)

**Subworkflows:**
- `subworkflows/pon_build.nf`: Implements PoN creation workflow
- `subworkflows/somatic_cnv.nf`: Placeholder for tumor sample processing (needs implementation)

**Configuration:**
- `nextflow.config`: Main configuration with profiles and parameters
- `conf/base.config`: Process resource configurations and publishing settings
- `conf/test.config`: Test-specific configurations

### Data Flow

1. **Sample parsing**: CSV samplesheet parsed to separate normal/tumor samples
2. **Interval processing**: Handles both BED files (preprocessed/annotated) and interval_list files
3. **PoN creation**: Normal samples → CollectReadCounts → CreateReadCountPanelOfNormals
4. **Tumor processing**: Uses PoN for denoising and segmentation (implementation pending)

### Input Requirements

**Samplesheet CSV columns:**
- `sample_id`: Unique sample identifier
- `type`: "normal" or "tumor"
- `cram`: Path to CRAM file
- `crai`: Path to CRAM index
- `sex`: Sample sex (optional)
- `tumor_normal_id`: Matched normal ID for tumor samples

**Required files:**
- Reference genome (FASTA, FAI, DICT)
- Intervals (BED or interval_list)
- SNP VCF for allelic counts (optional)

### Output Structure
- `results/pon/pon.hdf5`: Panel of Normals
- `results/denoise/<sample>/<sample>.denoisedCR.tsv`: Denoised copy ratios
- `results/segments/<sample>/<sample>.cr.seg`: Copy ratio segments
- `results/calls/<sample>/<sample>.called.seg`: Called segments

## Development Notes

- The `somatic_cnv.nf` subworkflow has been simplified to a basic placeholder structure to avoid compilation issues
- All GATK4 modules use nf-core standards except for local custom modules
- Resource configurations in `conf/base.config` are optimized for the workflow steps
- The workflow supports both WES and WGS data with appropriate parameter adjustments

## Recent Fixes Applied

**GATK Configuration Issues:**
- Added `ext.args = '--format HDF5 --interval-set-rule INTERSECTION'` to GATK4_COLLECTREADCOUNTS
- Added `ext.args = '--interval-merging-rule OVERLAPPING_ONLY'` to GATK4_PREPROCESSINTERVALS
- Changed `errorStrategy` from problematic closure to simple `'retry'`

**Workflow Structure Fixes:**
- Fixed channel handling in `pon_build.nf` to use proper `combine()` operations
- Fixed PREPROCESS_INTERVALS call structure to handle BED files properly
- Fixed tuple syntax in channel definitions (used `tuple()` wrapper)
- Fixed intervals channel extraction: `intervals_ch = preprocessed.interval_list.map { meta, interval_list -> interval_list }`

**Known Working State:**
- PoN building workflow is functional with BED format intervals
- PREPROCESS_INTERVALS successfully converts BED to interval_list format
- All channel operations and module calls have been debugged
- Workflow compiles and begins execution successfully

## Current Issues to Monitor

- Samplesheet paths must be accurate (fixed typos in sample paths)
- BED files work better than interval_list files with headers for this GATK version
- The workflow is currently set up for PoN-only mode (`--build_pon_only true`)

## File Format Requirements

**Samplesheet CSV:**
```csv
sample_id,type,cram,crai,sex,tumor_normal_id
ALK020-B-02,normal,/full/path/to/sample.recal.cram,/full/path/to/sample.recal.cram.crai,F
```

**Intervals:** 
- BED format works reliably (gets preprocessed automatically)
- interval_list format with headers can cause parsing issues in some GATK versions