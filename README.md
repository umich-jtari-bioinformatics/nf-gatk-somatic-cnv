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
├─ nf-tests/
│  ├─ nf-test.yml
│  ├─ subworkflows/
│  │  ├─ pon_build.test.nf
│  │  └─ somatic_cnv.test.nf
│  └─ data/
│     ├─ normals.csv
│     ├─ cohort.csv
│     ├─ tiny_intervals.bed
│     └─ tiny_snps.vcf.gz.tbi
└─ README.md
```


## Quick start

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
├── SOMATIC_CNV:GATK4_COLLECTALLELICCOUNTS_TUMOR
│   ├── ALK002.P.01.allelicCounts.tsv 
    └── ALK019.T.01.allelicCounts.tsv 

### Notes
	•	CRAMs must match the provided reference (including .dict & .fai)
	•	Intervals must be consistent across all steps (PoN and tumor)
	•	SNP VCF must match reference contigs; use common SNPs (gnomAD/1000G)