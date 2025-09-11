# nf-gatk-somatic-cnv
A Nextflow DSL2 workflow implementing the GATK somatic CNV pipeline:

`CollectReadCounts → CreateReadCountPanelOfNormals (PoN) → DenoiseReadCounts → CollectAllelicCounts → ModelSegments → CallCopyRatioSegments`

## Project Structure
```
nf-gatk-somatic-cnv/
├─ main.nf
├─ nextflow.config
├─ conf/
│  ├─ test.config
│  └─ base.config
├─ modules/
│  └─ gatk4/
│     ├─ collectreadcounts.nf
│     ├─ createreadcountpanelofnormals.nf
│     ├─ denoisereadcounts.nf            
│     ├─ collectalleliccounts.nf         # local
│     ├─ modelsegments.nf
│     ├─ callcopyratiosegments.nf        # local
│     ├─ preprocessintervals.nf          
│     └─ annotateintervals.nf            
├─ subworkflows/
│  ├─ pon_build.nf
│  └─ somatic_cnv.nf
├─ assets/
│  └─ README.assets.md
├─ schemas/
│  ├─ samplesheet.schema.json
│  └─ params.schema.json
├─ samples/
│  ├─ samplesheet_normals.csv
│  ├─ samplesheet_cohort.csv
│  ├─ tiny_ref.fa                        # stub (placeholder)
│  ├─ tiny_ref.fa.fai                    # stub (placeholder)
│  ├─ tiny_ref.dict                      # stub (placeholder)
│  ├─ tiny_intervals.bed                 # stub (placeholder)
│  └─ tiny_snps.vcf.gz.tbi               # stub (placeholder)
│
└─ README.md
```


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


## Quickstart

To make a PoN:
```bash
nextflow run /path/to/nf-gatk-somatic-cnv/main.nf \
  -profile docker \
  --samplesheet /home/ubuntu/code/nf-gatk-somatic-cnv/samples/samplesheet_normals.csv \
  --reference_fasta /home/ubuntu/code/nf-gatk-somatic-cnv/samples/Homo_sapiens_assembly38.fasta \
  --reference_fai /home/ubuntu/code/nf-gatk-somatic-cnv/samples/Homo_sapiens_assembly38.fasta.fai \
  --reference_dict /home/ubuntu/code/nf-gatk-somatic-cnv/samples/Homo_sapiens_assembly38.dict \
  --intervals /home/ubuntu/code/nf-gatk-somatic-cnv/samples/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.intervals \
  --snp_vcf /home/ubuntu/code/nf-gatk-somatic-cnv/samples/gnomad.snps.biallelic.norm.common.vcf.gz \
  --outdir ./results_pon \
  --build_pon_only true
```

-- For WES without intervals: add --assay wes --capture_bed /path/capture.bed (omit --intervals)

-- For WGS without intervals: add --assay wgs --bin_length 1000 --interval_padding 0 (omit --intervals)

To run a full tumor-only workflow using the PoN:

```bash
# assuming run from a project folder containing a config/ subfolder with a samplesheet
nextflow run /path/to/nf-gatk-somatic-cnv/main.nf \
   -profile singularity \
   --samplesheet ./config/samplesheet.gatk_cnv.non-alk.csv  \
   --reference_fasta /nfs/turbo/umms-jtari/jtsi-references/GATK/hg38/Homo_sapiens_assembly38.fasta \
   --reference_fai /nfs/turbo/umms-jtari/jtsi-references/GATK/hg38/Homo_sapiens_assembly38.fasta.fai \
   --reference_dict /nfs/turbo/umms-jtari/jtsi-references/GATK/hg38/Homo_sapiens_assembly38.dict \
   --intervals /nfs/turbo/umms-jtari/jtsi-references/GATK/somatic_cnv/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed \
   --snp_vcf /nfs/turbo/umms-jtari/jtsi-references/GATK/somatic_cnv/gnomad.snps.biallelic.norm.common.vcf.gz \
   --outdir ./results \
   --pon_hdf5 /nfs/turbo/umms-jtari/jtsi-references/GATK/somatic_cnv/pon.hdf5
```


### Inputs

See `samples/samplesheet_normals.csv` for required columns of the samplesheet.
Assumes indexed cram files of aligned reads
Requires matching reference.fasta reference.fasta.fai, and reference.dict 
Takes an intervals file in bed format (for WES data)
Needs a source of SNP loci, e.g. gnomad

### Outputs

#### For the PoN run
`results/pon/pon.hdf5`


#### Tumor samples
```
results/
├── denoise/
│   ├── ALK002.P.01_denoisedCR.tsv
│   └── ALK019.T.01_denoisedCR.tsv
├── read_counts
│   ├── ALK002.P.01.hdf5
│   └── ALK019.T.01.hdf5
├── segments
│   ├── ALK002.P.01.modelFinal.seg
│   ├── ALK019.T.01.modelFinal.seg
├── SOMATIC_CNV:GATK4_CALLCOPYRATIOSEGMENTS
│   ├── ALK002.P.01.calls.cns.tsv 
│   └── ALK019.T.01.calls.cns.tsv 
└── SOMATIC_CNV:GATK4_COLLECTALLELICCOUNTS_TUMOR
    ├── ALK002.P.01.allelicCounts.tsv 
    └── ALK019.T.01.allelicCounts.tsv 

```

### Notes
	•	CRAMs must match the provided reference (including .dict & .fai)
	•	Intervals must be consistent across all steps (PoN and tumor)
	•	SNP VCF must match reference contigs; use common SNPs (gnomAD/1000G)