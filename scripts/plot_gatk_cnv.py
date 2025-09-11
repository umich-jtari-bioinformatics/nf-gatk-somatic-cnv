#!/usr/bin/env python3
"""
Generate Sequenza-style genome plots from GATK CNV results.
Creates two-panel plots: BAF (top) and copy ratio (bottom).
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import argparse
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

def calculate_baf(ref_count, alt_count):
    """Calculate B-allele frequency from ref/alt counts"""
    total = ref_count + alt_count
    # Avoid division by zero
    baf = np.where(total > 0, alt_count / total, np.nan)
    return baf

def get_chromosome_positions(df, chr_lengths):
    """Convert chromosome coordinates to genome-wide positions"""
    df = df.copy()
    df['genome_start'] = 0
    df['genome_end'] = 0
    
    cumulative_pos = 0
    chr_boundaries = {}
    
    for chrom in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
                  'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
                  'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']:
        
        chr_data = df[df['CONTIG'] == chrom]
        if len(chr_data) == 0:
            continue
            
        chr_boundaries[chrom] = cumulative_pos
        
        # Update positions for this chromosome
        if 'START' in df.columns:
            df.loc[df['CONTIG'] == chrom, 'genome_start'] = df.loc[df['CONTIG'] == chrom, 'START'] + cumulative_pos
            df.loc[df['CONTIG'] == chrom, 'genome_end'] = df.loc[df['CONTIG'] == chrom, 'END'] + cumulative_pos
        else:
            df.loc[df['CONTIG'] == chrom, 'genome_start'] = df.loc[df['CONTIG'] == chrom, 'POSITION'] + cumulative_pos
            df.loc[df['CONTIG'] == chrom, 'genome_end'] = df.loc[df['CONTIG'] == chrom, 'POSITION'] + cumulative_pos
        
        # Move to next chromosome
        chr_length = chr_lengths.get(chrom, 250000000)  # Default length
        cumulative_pos += chr_length
    
    return df, chr_boundaries

def plot_gatk_cnv(copy_ratio_file, allelic_counts_file, segments_file, output_file, sample_id, model_segments_file=None):
    """Generate Sequenza-style CNV plot"""
    
    # Read data files
    print(f"Reading copy ratio data from {copy_ratio_file}")
    cr_data = read_gatk_file(copy_ratio_file)
    
    print(f"Reading allelic counts data from {allelic_counts_file}")
    ac_data = read_gatk_file(allelic_counts_file)
    
    print(f"Reading segments data from {segments_file}")
    seg_data = read_gatk_file(segments_file)
    
    # Read model segments for BAF segmentation if provided
    model_seg_data = None
    if model_segments_file and Path(model_segments_file).exists():
        print(f"Reading model segments data from {model_segments_file}")
        model_seg_data = read_gatk_file(model_segments_file)
    
    # Calculate BAF
    ac_data['BAF'] = calculate_baf(ac_data['REF_COUNT'], ac_data['ALT_COUNT'])
    # Filter out positions with no coverage
    ac_data = ac_data[(ac_data['REF_COUNT'] + ac_data['ALT_COUNT']) > 0]
    
    # Chromosome lengths (hg38 approximate)
    chr_lengths = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
    }
    
    # Convert to genome-wide coordinates
    cr_data, chr_boundaries = get_chromosome_positions(cr_data, chr_lengths)
    ac_data, _ = get_chromosome_positions(ac_data, chr_lengths)
    seg_data, _ = get_chromosome_positions(seg_data, chr_lengths)
    if model_seg_data is not None:
        model_seg_data, _ = get_chromosome_positions(model_seg_data, chr_lengths)
    
    # Create plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
    fig.suptitle(f'{sample_id} - GATK CNV Results', fontsize=14, fontweight='bold')
    
    # Top panel: BAF
    ax1.scatter(ac_data['genome_start'], ac_data['BAF'], 
               s=0.5, alpha=0.6, color='black')
    
    # Plot segmented BAF if available
    if model_seg_data is not None:
        for _, seg in model_seg_data.iterrows():
            if 'MINOR_ALLELE_FRACTION_POSTERIOR_50' in seg and not pd.isna(seg['MINOR_ALLELE_FRACTION_POSTERIOR_50']):
                baf_median = seg['MINOR_ALLELE_FRACTION_POSTERIOR_50']
                # Plot both minor allele fraction and 1-minor_allele_fraction (major allele)
                ax1.plot([seg['genome_start'], seg['genome_end']], [baf_median, baf_median], 
                        color='red', linewidth=2, alpha=0.8)
                ax1.plot([seg['genome_start'], seg['genome_end']], [1-baf_median, 1-baf_median], 
                        color='red', linewidth=2, alpha=0.8)
    
    ax1.set_ylabel('Allele frequency', fontsize=12)
    ax1.set_ylim(0, 1)
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0.5, color='red', linestyle='-', alpha=0.7, linewidth=1)
    
    # Bottom panel: Copy ratio
    # Plot individual points
    ax2.scatter(cr_data['genome_start'], 2**cr_data['LOG2_COPY_RATIO'], 
               s=0.5, alpha=0.6, color='black')
    
    # Plot segments
    for _, seg in seg_data.iterrows():
        if 'MEAN_LOG2_COPY_RATIO' in seg:
            ratio = 2**seg['MEAN_LOG2_COPY_RATIO']
            ax2.plot([seg['genome_start'], seg['genome_end']], [ratio, ratio], 
                    color='red', linewidth=2, alpha=0.8)
    
    ax2.set_ylabel('Depth ratio', fontsize=12)
    ax2.set_xlabel('Chromosome', fontsize=12)
    ax2.set_ylim(0, 2.5)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=1.0, color='red', linestyle='-', alpha=0.7, linewidth=1)
    
    # Add copy number annotations on right y-axis
    ax2_right = ax2.twinx()
    ax2_right.set_ylabel('Copy number', fontsize=12)
    ax2_right.set_ylim(0, 5)  # 0 to 5 copies
    ax2_right.set_yticks([0, 1, 2, 3, 4, 5])
    
    # Add chromosome boundaries and labels
    chr_centers = {}
    prev_pos = 0
    
    for chrom, start_pos in chr_boundaries.items():
        if chrom in chr_lengths:
            end_pos = start_pos + chr_lengths[chrom]
            center_pos = (start_pos + end_pos) / 2
            chr_centers[chrom] = center_pos
            
            # Add vertical lines at chromosome boundaries
            if start_pos > 0:
                ax1.axvline(x=start_pos, color='gray', alpha=0.5, linewidth=0.5)
                ax2.axvline(x=start_pos, color='gray', alpha=0.5, linewidth=0.5)
            
            prev_pos = end_pos
    
    # Set x-axis labels to chromosome names
    chr_labels = list(chr_centers.keys())
    chr_positions = list(chr_centers.values())
    
    # Only show labels for chromosomes that have data
    existing_chrs = []
    existing_pos = []
    for chr_name, pos in zip(chr_labels, chr_positions):
        if chr_name in cr_data['CONTIG'].unique():
            existing_chrs.append(chr_name.replace('chr', ''))
            existing_pos.append(pos)
    
    ax2.set_xticks(existing_pos)
    ax2.set_xticklabels(existing_chrs, rotation=0)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate Sequenza-style genome plots from GATK CNV results')
    parser.add_argument('--copy-ratio', required=True, help='Denoised copy ratio TSV file')
    parser.add_argument('--allelic-counts', required=True, help='Allelic counts TSV file')
    parser.add_argument('--segments', required=True, help='Copy ratio segments TSV file')
    parser.add_argument('--model-segments', help='Model segments TSV file (for segmented BAF)')
    parser.add_argument('--output', required=True, help='Output PNG file')
    parser.add_argument('--sample', required=True, help='Sample ID for plot title')
    
    args = parser.parse_args()
    
    # Check if files exist
    for filepath in [args.copy_ratio, args.allelic_counts, args.segments]:
        if not Path(filepath).exists():
            print(f"Error: File not found: {filepath}")
            sys.exit(1)
    
    plot_gatk_cnv(args.copy_ratio, args.allelic_counts, args.segments, args.output, args.sample, args.model_segments)

if __name__ == '__main__':
    main()