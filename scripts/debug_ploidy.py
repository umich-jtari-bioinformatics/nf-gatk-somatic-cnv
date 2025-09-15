#!/usr/bin/env python3
"""
Debug script to understand ploidy estimation issues
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

def read_gatk_file(filepath):
    """Read GATK TSV file, skipping SAM headers"""
    with open(filepath, 'r') as f:
        lines = []
        for line in f:
            if not line.startswith('@'):
                lines.append(line)
    
    from io import StringIO
    return pd.read_csv(StringIO(''.join(lines)), sep='\t')

def analyze_copy_ratios(segments_file, sample_id):
    """Analyze copy ratio distribution for a sample"""
    print(f"Analyzing {sample_id}")
    print("=" * 40)
    
    # Read segments
    segments_df = read_gatk_file(segments_file)
    
    # Filter to autosomal chromosomes
    autosomes = [f'chr{i}' for i in range(1, 23)]
    segments_auto = segments_df[segments_df['CONTIG'].isin(autosomes)]
    
    if 'MEAN_LOG2_COPY_RATIO' not in segments_auto.columns:
        print("No MEAN_LOG2_COPY_RATIO column found!")
        return
    
    # Calculate copy ratios and lengths
    segments_auto = segments_auto.copy()
    segments_auto['length'] = segments_auto['END'] - segments_auto['START']
    segments_auto['copy_ratio'] = 2 ** segments_auto['MEAN_LOG2_COPY_RATIO']
    
    # Show basic stats
    print(f"Total segments: {len(segments_auto)}")
    print(f"Total length: {segments_auto['length'].sum() / 1e6:.1f} Mb")
    print()
    
    # Length-weighted analysis
    total_length = segments_auto['length'].sum()
    segments_auto['weight'] = segments_auto['length'] / total_length
    
    weighted_mean_cr = np.sum(segments_auto['copy_ratio'] * segments_auto['weight'])
    median_cr = np.median(segments_auto['copy_ratio'])
    
    print(f"Weighted mean copy ratio: {weighted_mean_cr:.3f}")
    print(f"Median copy ratio: {median_cr:.3f}")
    print()
    
    # Create histogram
    copy_ratios = segments_auto['copy_ratio'].values
    weights = segments_auto['weight'].values
    
    # Show percentiles
    percentiles = [10, 25, 50, 75, 90]
    cr_percentiles = np.percentile(copy_ratios, percentiles)
    print("Copy ratio percentiles:")
    for p, cr in zip(percentiles, cr_percentiles):
        print(f"  {p}th: {cr:.3f}")
    print()
    
    # Find the mode using histogram
    hist, bin_edges = np.histogram(copy_ratios, bins=50, range=(0.1, 8.0), weights=weights)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    peak_idx = np.argmax(hist)
    mode_cr = bin_centers[peak_idx]
    
    print(f"Mode (peak) copy ratio: {mode_cr:.3f}")
    print(f"Current algorithm would estimate ploidy: {mode_cr * 2:.3f}")
    print()
    
    # Better ploidy estimation
    # The mode represents the most common copy number state
    # In a diploid normal, mode ≈ 1.0, so ploidy = 2
    # In a triploid tumor, mode ≈ 1.5, so ploidy = 3
    # So: estimated_ploidy = mode_cr * 2
    
    estimated_ploidy = mode_cr * 2
    print(f"Corrected ploidy estimate: {estimated_ploidy:.3f}")
    
    # Show distribution of calls
    if 'CALL' in segments_auto.columns:
        call_counts = segments_auto['CALL'].value_counts()
        call_lengths = segments_auto.groupby('CALL')['length'].sum() / 1e6
        print(f"\nCall distribution:")
        for call in ['-', '0', '+']:
            if call in call_counts.index:
                print(f"  {call}: {call_counts[call]} segments, {call_lengths[call]:.1f} Mb")
    
    return {
        'sample_id': sample_id,
        'weighted_mean_cr': weighted_mean_cr,
        'median_cr': median_cr,
        'mode_cr': mode_cr,
        'current_algorithm_ploidy': min(max(mode_cr * 2, 1.5), 6.0),
        'corrected_ploidy': mode_cr * 2
    }

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python debug_ploidy.py <segments_file> <sample_id>")
        sys.exit(1)
    
    segments_file = sys.argv[1]
    sample_id = sys.argv[2]
    
    if not Path(segments_file).exists():
        print(f"File not found: {segments_file}")
        sys.exit(1)
    
    result = analyze_copy_ratios(segments_file, sample_id)