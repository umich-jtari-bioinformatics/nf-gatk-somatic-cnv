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

echo "Starting gene CNV mapping for cohort..."
echo "Results directory: $RESULTS_DIR"
echo "Script: $SCRIPT_PATH"
echo "GENCODE GTF: $GENCODE_GTF"
echo "=========================================="

total_samples=0
processed_samples=0
failed_samples=0

mkdir -p gene_cnv_calls

for calls_file in "$RESULTS_DIR/calls_copyratiosegments"/*.calls.cns.tsv; do
    if [ ! -f "$calls_file" ]; then
        continue
    fi
    
    total_samples=$((total_samples + 1))
    
    basename_file=$(basename "$calls_file")
    sample_id="${basename_file%.calls.cns.tsv}"
    
    echo "Processing sample: $sample_id"
    
    output_file="gene_cnv_calls/${sample_id}_gene_cnv_calls.bed"
    
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
echo "Gene CNV mapping complete!"
echo "Total samples found: $total_samples"
echo "Successfully processed: $processed_samples"
echo "Failed: $failed_samples"
echo ""
echo "BED files saved to: gene_cnv_calls/"

if [ $failed_samples -gt 0 ]; then
    echo ""
    echo "Some samples failed processing. Check the output above for details."
    exit 1
fi

echo ""
echo "Creating summary statistics..."

summary_file="gene_cnv_calls/cohort_gene_cnv_summary.tsv"
echo -e "Sample_ID\tTotal_Genes\tNeutral_Genes\tAmplified_Genes\tDeleted_Genes\tNo_Coverage_Genes\tPercent_Amplified\tPercent_Deleted\tPercent_Altered" > "$summary_file"

for bed_file in gene_cnv_calls/*_gene_cnv_calls.bed; do
    if [ ! -f "$bed_file" ]; then
        continue
    fi
    
    sample_id=$(basename "$bed_file" _gene_cnv_calls.bed)
    
    # Count each call type (excluding header line)
    total_genes=$(grep -v "^track" "$bed_file" | wc -l | tr -d ' ')
    neutral_genes=$(grep -v "^track" "$bed_file" | grep ";neutral;" | wc -l | tr -d ' ')
    amplified_genes=$(grep -v "^track" "$bed_file" | grep ";amplified;" | wc -l | tr -d ' ')
    deleted_genes=$(grep -v "^track" "$bed_file" | grep ";deleted;" | wc -l | tr -d ' ')
    no_coverage_genes=$(grep -v "^track" "$bed_file" | grep ";no_coverage" | wc -l | tr -d ' ')
    
    # Calculate percentages
    if [ "$total_genes" -gt 0 ]; then
        percent_amplified=$(echo "scale=1; $amplified_genes * 100 / $total_genes" | bc)
        percent_deleted=$(echo "scale=1; $deleted_genes * 100 / $total_genes" | bc)
        percent_altered=$(echo "scale=1; ($amplified_genes + $deleted_genes) * 100 / $total_genes" | bc)
    else
        percent_amplified="0.0"
        percent_deleted="0.0" 
        percent_altered="0.0"
    fi
    
    echo -e "$sample_id\t$total_genes\t$neutral_genes\t$amplified_genes\t$deleted_genes\t$no_coverage_genes\t$percent_amplified\t$percent_deleted\t$percent_altered" >> "$summary_file"
done

echo "Summary statistics created: $summary_file"
echo ""
echo "Usage tips:"
echo "  - Load BED files into IGV or UCSC Genome Browser for visualization"
echo "  - Use summary file to identify samples with high gene alteration rates"
echo "  - Search BED files for specific genes of interest (e.g., ALK, EGFR, TP53)"
echo ""
echo "Example searches:"
echo "  grep 'ALK;' gene_cnv_calls/*.bed | head -10"
echo "  grep 'amplified' gene_cnv_calls/ALK019.T.01_gene_cnv_calls.bed | head -5"