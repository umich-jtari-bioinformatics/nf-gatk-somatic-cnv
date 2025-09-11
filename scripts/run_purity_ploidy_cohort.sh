#!/bin/bash

set -euo pipefail

RESULTS_DIR="/Users/pulintz/Data/JTSI/Sequencing/jtari_patient_sequencing/nf-gatk-somatic-cnv/results"
SCRIPT_DIR="~/Code/nf-gatk-somatic-cnv/scripts"
SCRIPT_PATH="$SCRIPT_DIR/estimate_purity_ploidy.py"

if [ ! -d "$RESULTS_DIR" ]; then
    echo "Error: Results directory not found: $RESULTS_DIR"
    exit 1
fi

if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: Script not found: $SCRIPT_PATH"
    exit 1
fi

echo "Starting purity/ploidy analysis for cohort..."
echo "Results directory: $RESULTS_DIR"
echo "Script: $SCRIPT_PATH"
echo "=========================================="

total_samples=0
processed_samples=0
failed_samples=0

mkdir -p purity_ploidy_reports

for denoise_file in "$RESULTS_DIR/denoise"/*.tsv; do
    if [ ! -f "$denoise_file" ]; then
        continue
    fi
    
    total_samples=$((total_samples + 1))
    
    basename_file=$(basename "$denoise_file")
    sample_id="${basename_file%_denoisedCR.tsv}"
    
    echo "Processing sample: $sample_id"
    
    copy_ratio_file="$denoise_file"
    segments_file="$RESULTS_DIR/calls_copyratiosegments/${sample_id}.calls.cns.tsv"
    model_segments_file="$RESULTS_DIR/segments/${sample_id}.modelFinal.seg"
    output_file="purity_ploidy_reports/${sample_id}_purity_ploidy_report.txt"
    
    if [ ! -f "$copy_ratio_file" ]; then
        echo "  Warning: Copy ratio file not found: $copy_ratio_file"
    fi
    
    if [ ! -f "$segments_file" ]; then
        echo "  Warning: Segments file not found: $segments_file"
        segments_file=""
    fi
    
    if [ ! -f "$model_segments_file" ]; then
        echo "  Warning: Model segments file not found: $model_segments_file"
        model_segments_file=""
    fi
    
    cmd_args="--sample $sample_id --output $output_file"
    
    if [ -f "$copy_ratio_file" ]; then
        cmd_args="$cmd_args --copy-ratio $copy_ratio_file"
    fi
    
    if [ -n "$segments_file" ] && [ -f "$segments_file" ]; then
        cmd_args="$cmd_args --segments $segments_file"
    fi
    
    if [ -n "$model_segments_file" ] && [ -f "$model_segments_file" ]; then
        cmd_args="$cmd_args --model-segments $model_segments_file"
    fi
    
    if python3 "$SCRIPT_PATH" $cmd_args; then
        processed_samples=$((processed_samples + 1))
        echo "  ✓ Successfully processed $sample_id"
    else
        failed_samples=$((failed_samples + 1))
        echo "  ✗ Failed to process $sample_id"
    fi
    
    echo ""
done

echo "=========================================="
echo "Analysis complete!"
echo "Total samples found: $total_samples"
echo "Successfully processed: $processed_samples"
echo "Failed: $failed_samples"
echo ""
echo "Reports saved to: purity_ploidy_reports/"

if [ $failed_samples -gt 0 ]; then
    echo ""
    echo "Some samples failed processing. Check the output above for details."
    exit 1
fi

echo ""
echo "Creating summary report..."

summary_file="purity_ploidy_reports/cohort_summary.tsv"
echo -e "Sample_ID\tEstimated_Purity\tPurity_Range_Min\tPurity_Range_Max\tEstimated_Ploidy\tConfidence\tMedian_Copy_Ratio" > "$summary_file"

for report_file in purity_ploidy_reports/*_purity_ploidy_report.txt; do
    if [ ! -f "$report_file" ]; then
        continue
    fi
    
    sample_id=$(basename "$report_file" _purity_ploidy_report.txt)
    
    purity=$(grep "Estimated tumor purity:" "$report_file" | cut -d: -f2 | tr -d ' ' || echo "NA")
    purity_range=$(grep "Purity range:" "$report_file" | cut -d: -f2 | tr -d ' ' || echo "NA")
    ploidy=$(grep "Estimated tumor ploidy:" "$report_file" | cut -d: -f2 | tr -d ' ' || echo "NA")
    confidence=$(grep "Confidence:" "$report_file" | cut -d: -f2 | tr -d ' ' || echo "NA")
    median_cr=$(grep "Median copy ratio:" "$report_file" | cut -d: -f2 | tr -d ' ' || echo "NA")
    
    if [ "$purity_range" != "NA" ] && [[ "$purity_range" == *" - "* ]]; then
        purity_min=$(echo "$purity_range" | cut -d' ' -f1)
        purity_max=$(echo "$purity_range" | cut -d' ' -f3)
    else
        purity_min="NA"
        purity_max="NA"
    fi
    
    echo -e "$sample_id\t$purity\t$purity_min\t$purity_max\t$ploidy\t$confidence\t$median_cr" >> "$summary_file"
done

echo "Summary report created: $summary_file"
echo ""
echo "To review results:"
echo "  - Individual reports: purity_ploidy_reports/*_purity_ploidy_report.txt"
echo "  - Summary table: $summary_file"