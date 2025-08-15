# nf-gatk-somatic-cnv
Implementation of the GATK somatic CNV workflow for WES and WGS data

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



