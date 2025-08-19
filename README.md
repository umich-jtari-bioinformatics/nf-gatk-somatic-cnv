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
│  ├─ tiny_ref.fa                        # stub (placeholder text)
│  ├─ tiny_ref.fa.fai                    # stub (placeholder text)
│  ├─ tiny_ref.dict                      # stub (placeholder text)
│  ├─ tiny_intervals.bed                 # stub
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

```bash
nextflow run main.nf -profile docker \
  --samplesheet assets/samplesheet_full.csv \
  --reference_fasta /path/ref.fa --reference_fai /path/ref.fa.fai --reference_dict /path/ref.dict \
  --intervals /path/intervals.interval_list \
  --snp_vcf /path/common_snps.vcf.gz \
  --outdir results
```

-- For WES without intervals: add --assay wes --capture_bed /path/capture.bed (omit --intervals)

-- For WGS without intervals: add --assay wgs --bin_length 1000 --interval_padding 0 (omit --intervals)

### Inputs

See `assets/samplesheet_full.csv` and `schemas/samplesheet_schema.json` for required columns.

### Outputs
`results/pon/pon.hdf5`

Per tumor:

-- `results/denoise/<sample>/<sample>.denoisedCR.tsv`

-- `results/segments/<sample>/<sample>.cr.seg`

-- `results/calls/<sample>/<sample>.called.seg`

### Notes
	•	CRAMs must match the provided reference (including .dict & .fai)
	•	Intervals must be consistent across all steps (PoN and tumor)
	•	SNP VCF must match reference contigs; use common SNPs (gnomAD/1000G)