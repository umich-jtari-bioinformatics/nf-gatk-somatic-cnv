#!/bin/bash

set -euo pipefail

RESULTS_DIR="/Users/pulintz/Data/JTSI/Sequencing/jtari_patient_sequencing/nf-gatk-somatic-cnv/results"
SCRIPT_DIR="/Users/pulintz/Code/nf-gatk-somatic-cnv/scripts"
GENCODE_GTF="/Users/pulintz/Data/references/hg38/Gencode/gencode.v49.annotation.mainChr.gtf.gz"
SCRIPT_PATH="$SCRIPT_DIR/genes_to_cnv_calls.py"

if [ ! -d "$RESULTS_DIR" ]; then
    echo "Error: Results directory not found: $RESULTS_DIR"
    exit 1
fi

if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: Script not found: $SCRIPT_PATH"
    exit 1
fi

if [ ! -f "$GENCODE_GTF" ]; then
    echo "Error: GENCODE GTF file not found: $GENCODE_GTF"
    exit 1
fi

echo "Starting gene CNV VCF generation for cohort..."
echo "Results directory: $RESULTS_DIR"
echo "Script: $SCRIPT_PATH"
echo "GENCODE GTF: $GENCODE_GTF"
echo "=========================================="

total_samples=0
processed_samples=0
failed_samples=0

mkdir -p gene_cnv_vcfs

for calls_file in "$RESULTS_DIR/calls_copyratiosegments"/*.calls.cns.tsv; do
    if [ ! -f "$calls_file" ]; then
        continue
    fi
    
    total_samples=$((total_samples + 1))
    
    basename_file=$(basename "$calls_file")
    sample_id="${basename_file%.calls.cns.tsv}"
    
    echo "Processing sample: $sample_id"
    
    output_file="gene_cnv_vcfs/${sample_id}_gene_cnv_calls.vcf"
    
    if python3 "$SCRIPT_PATH" \
        --gtf "$GENCODE_GTF" \
        --calls "$calls_file" \
        --sample "$sample_id" \
        --output "$output_file"; then
        processed_samples=$((processed_samples + 1))
        echo "  ✓ Successfully processed $sample_id"
    else
        failed_samples=$((failed_samples + 1))
        echo "  ✗ Failed to process $sample_id"
    fi
    
    echo ""
done

echo "=========================================="
echo "Gene CNV VCF generation complete!"
echo "Total samples found: $total_samples"
echo "Successfully processed: $processed_samples"
echo "Failed: $failed_samples"
echo ""
echo "VCF files saved to: gene_cnv_vcfs/"

if [ $failed_samples -gt 0 ]; then
    echo ""
    echo "Some samples failed processing. Check the output above for details."
    exit 1
fi

echo ""
echo "Creating summary statistics..."

summary_file="gene_cnv_vcfs/cohort_gene_cnv_summary.tsv"
echo -e "Sample_ID\tTotal_CNV_Calls\tAmplifications\tDeletions\tPercent_Amplified\tPercent_Deleted" > "$summary_file"

for vcf_file in gene_cnv_vcfs/*_gene_cnv_calls.vcf; do
    if [ ! -f "$vcf_file" ]; then
        continue
    fi
    
    sample_id=$(basename "$vcf_file" _gene_cnv_calls.vcf)
    
    # Count CNV calls (excluding header lines)
    total_cnvs=$(grep -v "^#" "$vcf_file" | wc -l | tr -d ' ')
    amplifications=$(grep -v "^#" "$vcf_file" | grep "SVTYPE=DUP" | wc -l | tr -d ' ')
    deletions=$(grep -v "^#" "$vcf_file" | grep "SVTYPE=DEL" | wc -l | tr -d ' ')
    
    # Calculate percentages
    if [ "$total_cnvs" -gt 0 ]; then
        percent_amplified=$(echo "scale=1; $amplifications * 100 / $total_cnvs" | bc)
        percent_deleted=$(echo "scale=1; $deletions * 100 / $total_cnvs" | bc)
    else
        percent_amplified="0.0"
        percent_deleted="0.0"
    fi
    
    echo -e "$sample_id\t$total_cnvs\t$amplifications\t$deletions\t$percent_amplified\t$percent_deleted" >> "$summary_file"
done

echo "Summary statistics created: $summary_file"
echo ""
echo "Usage tips:"
echo "  - Use bcftools/vcftools for VCF manipulation and analysis"
echo "  - Load VCF files into IGV for visualization"
echo "  - Query specific genes: bcftools query -f '%CHROM\\t%POS\\t%INFO/GENE\\t%INFO/CNV_CALL\\n' file.vcf"
echo ""
echo "Example VCF queries:"
echo "  bcftools view -H gene_cnv_vcfs/ALK019.T.01_gene_cnv_calls.vcf | grep 'GENE=ALK'"
echo "  bcftools query -f '%INFO/GENE\\t%INFO/CNV_CALL\\t%INFO/LOG2CR\\n' gene_cnv_vcfs/*.vcf | grep -E 'ALK|EGFR|TP53'"