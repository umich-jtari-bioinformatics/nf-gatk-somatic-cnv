#!/usr/bin/env python3
"""
Estimate tumor purity and ploidy from GATK CNV results.

This script analyzes copy ratio and allelic fraction data to estimate:
1. Tumor purity (cellularity)
2. Average tumor ploidy
3. Confidence metrics for the estimates

Based on the principle that:
- Pure diploid samples have copy ratio centered at 1.0 (log2 = 0)
- Contamination with normal cells shifts copy ratios toward 1.0
- BAF patterns in heterozygous deletions reveal purity
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path
from collections import Counter

def read_gatk_file(filepath):
    """Read GATK TSV file, skipping SAM headers"""
    with open(filepath, 'r') as f:
        lines = []
        for line in f:
            if not line.startswith('@'):
                lines.append(line)
    
    from io import StringIO
    return pd.read_csv(StringIO(''.join(lines)), sep='\t')

def filter_autosomal_chromosomes(df):
    """Filter to autosomal chromosomes only"""
    autosomes = [f'chr{i}' for i in range(1, 23)]
    return df[df['CONTIG'].isin(autosomes)]

def estimate_baseline_ploidy(segments_df, copy_ratio_df):
    """
    Estimate baseline tumor ploidy accounting for GATK PoN normalization.
    
    Key insight: GATK normalizes to diploid PoN, so copy ratios represent
    deviations from tumor baseline, not absolute copy numbers.
    Use amplification/deletion balance to infer true baseline ploidy.
    """
    # Filter to autosomal chromosomes
    segments_auto = filter_autosomal_chromosomes(segments_df)
    
    if 'MEAN_LOG2_COPY_RATIO' in segments_auto.columns and 'CALL' in segments_auto.columns:
        # Weight by segment length
        segments_auto = segments_auto.copy()
        segments_auto['length'] = segments_auto['END'] - segments_auto['START']
        segments_auto['copy_ratio'] = 2 ** segments_auto['MEAN_LOG2_COPY_RATIO']
        
        total_length = segments_auto['length'].sum()
        
        # Analyze call distribution to infer baseline ploidy
        call_fractions = {}
        call_mean_cr = {}
        
        for call in ['-', '0', '+']:
            call_segs = segments_auto[segments_auto['CALL'] == call]
            if len(call_segs) > 0:
                call_fractions[call] = call_segs['length'].sum() / total_length
                call_mean_cr[call] = call_segs['copy_ratio'].mean()
            else:
                call_fractions[call] = 0.0
                call_mean_cr[call] = 1.0
        
        amp_fraction = call_fractions['+']
        del_fraction = call_fractions['-']
        neutral_fraction = call_fractions['0']
        neutral_cr = call_mean_cr['0']
        
        # Method 1: Amplification/deletion balance
        # Strong excess of amplifications suggests higher baseline ploidy
        amp_del_ratio = amp_fraction / (del_fraction + 0.01)
        
        ploidy_estimates = []
        
        # Estimate from amp/del imbalance
        if amp_del_ratio > 1.5 and amp_fraction > 0.1:
            # Excess amplifications suggest baseline > 2N
            balance_ploidy = 2.0 + min(amp_del_ratio - 1.0, 3.0) * 0.5
            ploidy_estimates.append(balance_ploidy)
        
        # Method 2: Neutral region analysis
        # If neutral regions show copy ratio != 1.0, adjust baseline
        if neutral_fraction > 0.3:  # Sufficient neutral regions
            if neutral_cr < 0.95:
                # Neutral regions appear deleted relative to PoN
                # Suggests baseline > 2N (PoN thinks it's deletion)
                cr_ploidy = 2.0 + (1.0 - neutral_cr) * 2.0
                ploidy_estimates.append(cr_ploidy)
            elif neutral_cr > 1.05:
                # Neutral regions appear amplified relative to PoN
                # Suggests baseline < 2N (unlikely in cancer)
                cr_ploidy = 2.0 * neutral_cr
                ploidy_estimates.append(cr_ploidy)
        
        # Method 3: Traditional histogram approach (as fallback)
        weighted_ratios = []
        for _, seg in segments_auto.iterrows():
            weight = max(1, int(seg['length'] / 1000000))
            weighted_ratios.extend([seg['copy_ratio']] * weight)
        
        if weighted_ratios:
            ratios_array = np.array(weighted_ratios)
            hist, bin_edges = np.histogram(ratios_array, bins=50, range=(0.3, 6.0))
            peak_idx = np.argmax(hist)
            peak_copy_ratio = (bin_edges[peak_idx] + bin_edges[peak_idx + 1]) / 2
            
            # For GATK data, peak represents tumor baseline relative to diploid PoN
            hist_ploidy = 2.0 * peak_copy_ratio
            ploidy_estimates.append(hist_ploidy)
        
        # Consensus ploidy estimate
        if ploidy_estimates:
            baseline_ploidy = np.median(ploidy_estimates)
            # Apply reasonable bounds but allow higher ploidy
            return max(1.5, min(baseline_ploidy, 8.0))
        else:
            return 2.0
    
    return 2.0  # Default diploid

def estimate_purity_from_baf(model_segments_df, min_segment_length=2e6):
    """
    Estimate purity from BAF patterns in copy-neutral LOH and deletions.
    
    In pure tumors:
    - Heterozygous deletions show BAF ~0 or ~1
    - Copy-neutral LOH shows BAF ~0 or ~1
    - Heterozygous regions show BAF ~0.5
    
    In impure tumors, these values are shifted toward 0.5
    """
    if model_segments_df is None or model_segments_df.empty:
        return None, None
    
    segments_auto = filter_autosomal_chromosomes(model_segments_df)
    segments_auto = segments_auto.copy()
    segments_auto['length'] = segments_auto['END'] - segments_auto['START']
    
    # Filter segments by minimum length
    long_segments = segments_auto[segments_auto['length'] >= min_segment_length]
    
    if long_segments.empty:
        return None, None
    
    purity_estimates = []
    
    # Method 1: Look for copy-neutral LOH (copy ratio ~1, BAF ~0 or ~1)
    if 'LOG2_COPY_RATIO_POSTERIOR_50' in long_segments.columns and 'MINOR_ALLELE_FRACTION_POSTERIOR_50' in long_segments.columns:
        for _, seg in long_segments.iterrows():
            log2_cr = seg['LOG2_COPY_RATIO_POSTERIOR_50']
            maf = seg['MINOR_ALLELE_FRACTION_POSTERIOR_50']
            
            if pd.isna(log2_cr) or pd.isna(maf):
                continue
                
            copy_ratio = 2 ** log2_cr
            
            # Copy-neutral segments (copy ratio 0.7-1.3)
            if 0.7 <= copy_ratio <= 1.3:
                # If this segment shows LOH (low minor allele fraction)
                if maf < 0.45:  # Should be ~0.5 in heterozygous diploid
                    # In pure tumor: maf would be ~0
                    # In impure tumor: maf = purity * 0 + (1-purity) * 0.5 = (1-purity) * 0.5
                    # So: purity = 1 - 2*maf (when expected pure maf = 0)
                    purity_est = 1 - 2 * maf
                    if 0.3 <= purity_est <= 1.0:  # Reasonable range
                        purity_estimates.append(purity_est)
    
    # Method 2: Look at heterozygous deletions (copy ratio ~0.5, BAF should be 0 or 1 in pure)
    for _, seg in long_segments.iterrows():
        log2_cr = seg.get('LOG2_COPY_RATIO_POSTERIOR_50')
        maf = seg.get('MINOR_ALLELE_FRACTION_POSTERIOR_50')
        
        if pd.isna(log2_cr) or pd.isna(maf):
            continue
            
        copy_ratio = 2 ** log2_cr
        
        # Heterozygous deletion (copy ratio 0.3-0.8)
        if 0.3 <= copy_ratio <= 0.8:
            # In pure tumor, remaining allele should be either A or B (maf ~0 or ~1)
            # In impure tumor, gets diluted toward 0.5
            deviation_from_half = abs(maf - 0.5)
            if deviation_from_half > 0.05:  # Some signal
                # Estimate purity from how far from 0.5 the BAF is
                # Pure tumor: maf = 0 or 1
                # Impure: maf = purity * (0 or 1) + (1-purity) * 0.5
                expected_pure_maf = 0 if maf < 0.5 else 1
                # maf = purity * expected_pure_maf + (1-purity) * 0.5
                # Solving: purity = 2 * abs(maf - 0.5)
                purity_est = 2 * deviation_from_half
                if 0.3 <= purity_est <= 1.0:
                    purity_estimates.append(purity_est)
    
    if purity_estimates:
        purity_mean = np.mean(purity_estimates)
        purity_std = np.std(purity_estimates)
        return purity_mean, purity_std
    
    return None, None

def estimate_purity_from_copy_ratio(segments_df, estimated_ploidy=2.0):
    """
    Estimate purity from copy ratio shifts.
    
    In a pure tumor, copy ratios reflect true copy numbers.
    In impure tumors, copy ratios are shifted toward normal (diploid).
    """
    segments_auto = filter_autosomal_chromosomes(segments_df)
    
    if 'MEAN_LOG2_COPY_RATIO' not in segments_auto.columns:
        return None, None
    
    segments_auto = segments_auto.copy()
    segments_auto['length'] = segments_auto['END'] - segments_auto['START']
    segments_auto['copy_ratio'] = 2 ** segments_auto['MEAN_LOG2_COPY_RATIO']
    
    # Look for clear amplifications (copy ratio > 1.3) and deletions (copy ratio < 0.8)
    amplifications = segments_auto[segments_auto['copy_ratio'] > 1.3]
    deletions = segments_auto[segments_auto['copy_ratio'] < 0.8]
    
    purity_estimates = []
    
    # For amplifications: observed = purity * true + (1-purity) * 2
    for _, seg in amplifications.iterrows():
        if seg['length'] > 2e6:  # Minimum 2Mb
            observed_cr = seg['copy_ratio']
            # Assume true copy number in pure tumor would be 3 or 4
            for true_cn in [3, 4, 5, 6]:
                # observed = purity * true_cn + (1-purity) * 2
                # purity = (observed - 2) / (true_cn - 2)
                if true_cn > 2:
                    purity_est = (observed_cr - 2) / (true_cn - 2)
                    if 0.3 <= purity_est <= 1.0:
                        purity_estimates.append(purity_est)
    
    # For deletions: similar logic
    for _, seg in deletions.iterrows():
        if seg['length'] > 2e6:
            observed_cr = seg['copy_ratio']
            # Assume true copy number in pure tumor would be 1 or 0
            for true_cn in [0, 1]:
                # observed = purity * true_cn + (1-purity) * 2
                # purity = (observed - 2) / (true_cn - 2)
                if true_cn < 2:
                    purity_est = (2 - observed_cr) / (2 - true_cn)
                    if 0.3 <= purity_est <= 1.0:
                        purity_estimates.append(purity_est)
    
    if purity_estimates:
        purity_mean = np.mean(purity_estimates)
        purity_std = np.std(purity_estimates)
        return purity_mean, purity_std
    
    return None, None

def estimate_median_copy_ratio_shift(copy_ratio_df):
    """
    Estimate purity from overall shift in copy ratio distribution.
    Pure diploid should be centered at 1.0, impurity shifts toward 1.0.
    """
    if copy_ratio_df is None or copy_ratio_df.empty:
        return None
    
    auto_data = filter_autosomal_chromosomes(copy_ratio_df)
    
    if 'LOG2_COPY_RATIO' not in auto_data.columns:
        return None
    
    # Convert log2 ratios to linear ratios
    copy_ratios = 2 ** auto_data['LOG2_COPY_RATIO']
    
    # Remove extreme outliers
    q1, q99 = np.percentile(copy_ratios, [1, 99])
    filtered_ratios = copy_ratios[(copy_ratios >= q1) & (copy_ratios <= q99)]
    
    if len(filtered_ratios) < 100:
        return None
    
    median_ratio = np.median(filtered_ratios)
    return median_ratio

def generate_report(copy_ratio_file, segments_file, model_segments_file, sample_id, output_file):
    """Generate purity/ploidy estimation report"""
    
    print(f"Analyzing sample: {sample_id}")
    print("=" * 50)
    
    # Read data
    copy_ratio_df = None
    if copy_ratio_file and Path(copy_ratio_file).exists():
        print(f"Reading copy ratio data: {copy_ratio_file}")
        copy_ratio_df = read_gatk_file(copy_ratio_file)
    
    segments_df = None
    if segments_file and Path(segments_file).exists():
        print(f"Reading segments data: {segments_file}")
        segments_df = read_gatk_file(segments_file)
    
    model_segments_df = None
    if model_segments_file and Path(model_segments_file).exists():
        print(f"Reading model segments data: {model_segments_file}")
        model_segments_df = read_gatk_file(model_segments_file)
    
    print("\nEstimation Results:")
    print("-" * 30)
    
    # Estimate ploidy
    estimated_ploidy = 2.0
    if segments_df is not None:
        estimated_ploidy = estimate_baseline_ploidy(segments_df, copy_ratio_df)
        print(f"Estimated tumor ploidy: {estimated_ploidy:.2f}")
    else:
        print("Estimated tumor ploidy: 2.0 (default - no segments data)")
    
    # Estimate purity using multiple methods
    purity_estimates = []
    methods = []
    failed_methods = []
    
    # Method 1: BAF-based estimation
    if model_segments_df is not None:
        purity_baf, purity_baf_std = estimate_purity_from_baf(model_segments_df)
        if purity_baf is not None:
            purity_estimates.append(purity_baf)
            methods.append(f"BAF analysis")
            print(f"Purity from BAF analysis: {purity_baf:.2f} ± {purity_baf_std:.2f}")
        else:
            failed_methods.append("BAF: No suitable segments found")
            print("Purity from BAF analysis: Unable to estimate")
    else:
        failed_methods.append("BAF: Model segments file missing")
    
    # Method 2: Copy ratio shift
    if segments_df is not None:
        purity_cr, purity_cr_std = estimate_purity_from_copy_ratio(segments_df, estimated_ploidy)
        if purity_cr is not None:
            purity_estimates.append(purity_cr)
            methods.append("Copy ratio analysis")
            print(f"Purity from copy ratio analysis: {purity_cr:.2f} ± {purity_cr_std:.2f}")
        else:
            failed_methods.append("Copy ratio: No clear CNAs found")
            print("Purity from copy ratio analysis: Unable to estimate")
    else:
        failed_methods.append("Copy ratio: Segments file missing")
    
    # Method 3: Overall distribution shift
    median_cr = None
    if copy_ratio_df is not None:
        median_cr = estimate_median_copy_ratio_shift(copy_ratio_df)
        if median_cr is not None:
            print(f"Median copy ratio: {median_cr:.3f}")
            # Simple heuristic: deviation from 1.0 suggests purity issues
            if abs(median_cr - 1.0) < 0.1:
                print("Overall copy ratio distribution: Close to diploid baseline")
            else:
                print(f"Overall copy ratio distribution: Shifted from diploid baseline")
    
    # Calculate fraction of genome with copy number alterations
    genome_alteration_fraction = None
    if segments_df is not None:
        segments_auto = filter_autosomal_chromosomes(segments_df)
        if 'CALL' in segments_auto.columns:
            segments_auto = segments_auto.copy()
            segments_auto['length'] = segments_auto['END'] - segments_auto['START']
            total_length = segments_auto['length'].sum()
            
            # Calculate fractions for each call type
            call_fractions = {}
            for call in ['-', '0', '+']:
                call_segs = segments_auto[segments_auto['CALL'] == call]
                call_fractions[call] = call_segs['length'].sum() / total_length if len(call_segs) > 0 else 0
            
            genome_alteration_fraction = call_fractions['+'] + call_fractions['-']
            
            print(f"Genome composition: {call_fractions['+']:.1%} amplified, {call_fractions['0']:.1%} neutral, {call_fractions['-']:.1%} deleted")
            print(f"Total genome altered: {genome_alteration_fraction:.1%}")
            
            # Assess genomic instability level
            if genome_alteration_fraction > 0.4:
                instability_level = "High"
            elif genome_alteration_fraction > 0.15:
                instability_level = "Moderate"  
            elif genome_alteration_fraction > 0.05:
                instability_level = "Low"
            else:
                instability_level = "Minimal"
            print(f"Genomic instability: {instability_level}")
    
    # Method 4: Fallback purity estimation from CALL distribution
    # If no clear CNAs found by primary methods, look at overall genomic instability
    if not purity_estimates and segments_df is not None:
        segments_auto = filter_autosomal_chromosomes(segments_df)
        if 'CALL' in segments_auto.columns:
            segments_auto = segments_auto.copy()
            segments_auto['length'] = segments_auto['END'] - segments_auto['START']
            total_length = segments_auto['length'].sum()
            
            # Calculate fraction of genome with any CNA
            altered_segs = segments_auto[segments_auto['CALL'] != '0']
            altered_fraction = altered_segs['length'].sum() / total_length if len(altered_segs) > 0 else 0
            
            # If significant portion of genome shows CNAs, estimate purity from instability
            if altered_fraction > 0.05:  # At least 5% of genome altered
                # More instability suggests higher purity
                instability_purity = min(0.9, 0.3 + altered_fraction * 0.8)
                purity_estimates.append(instability_purity)
                methods.append("Genomic instability")
                print(f"Purity from genomic instability: {instability_purity:.2f} ({altered_fraction:.1%} genome altered)")
            else:
                failed_methods.append("Instability: Low genomic alteration")
                print(f"Purity from genomic instability: Unable to estimate ({altered_fraction:.1%} genome altered)")
    
    # Consensus estimate
    if purity_estimates:
        consensus_purity = np.median(purity_estimates)
        purity_range = [np.min(purity_estimates), np.max(purity_estimates)]
        print(f"\nConsensus purity estimate: {consensus_purity:.2f}")
        print(f"Range across methods: {purity_range[0]:.2f} - {purity_range[1]:.2f}")
        
        # Quality assessment
        purity_std = np.std(purity_estimates)
        if purity_std < 0.1:
            quality = "High"
        elif purity_std < 0.2:
            quality = "Medium"
        else:
            quality = "Low"
        print(f"Estimate confidence: {quality} (std dev: {purity_std:.2f})")
    else:
        consensus_purity = None
        print("\nUnable to estimate purity from available data")
        quality = "Unable to estimate"
    
    # Write report
    if output_file:
        with open(output_file, 'w') as f:
            f.write(f"Purity/Ploidy Analysis Report\n")
            f.write(f"Sample: {sample_id}\n")
            f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write(f"Estimated tumor ploidy: {estimated_ploidy:.2f}\n")
            
            if consensus_purity is not None:
                f.write(f"Estimated tumor purity: {consensus_purity:.2f}\n")
                f.write(f"Purity range: {purity_range[0]:.2f} - {purity_range[1]:.2f}\n")
                f.write(f"Confidence: {quality}\n")
            else:
                f.write(f"Estimated tumor purity: Unable to estimate\n")
            
            f.write(f"\nMethods used: {', '.join(methods) if methods else 'None'}\n")
            if failed_methods:
                f.write(f"Failed methods: {'; '.join(failed_methods)}\n")
            
            if median_cr is not None:
                f.write(f"Median copy ratio: {median_cr:.3f}\n")
            
            if genome_alteration_fraction is not None:
                f.write(f"Genome altered: {genome_alteration_fraction:.1%}\n")
                f.write(f"Genomic instability: {instability_level}\n")
        
        print(f"\nReport saved to: {output_file}")
    
    return {
        'sample_id': sample_id,
        'estimated_ploidy': estimated_ploidy,
        'estimated_purity': consensus_purity,
        'purity_range': purity_range if purity_estimates else None,
        'confidence': quality if purity_estimates else 'Unable to estimate',
        'median_copy_ratio': median_cr,
        'genome_alteration_fraction': genome_alteration_fraction,
        'genomic_instability': instability_level if genome_alteration_fraction is not None else None,
        'methods_used': ', '.join(methods) if methods else 'None',
        'failed_methods': '; '.join(failed_methods) if failed_methods else 'None'
    }

def main():
    parser = argparse.ArgumentParser(description='Estimate tumor purity and ploidy from GATK CNV results')
    parser.add_argument('--copy-ratio', help='Denoised copy ratio TSV file')
    parser.add_argument('--segments', help='Copy ratio segments TSV file (calls.cns.tsv)')
    parser.add_argument('--model-segments', help='Model segments TSV file (modelFinal.seg)')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--output', help='Output report file')
    
    args = parser.parse_args()
    
    if not any([args.copy_ratio, args.segments, args.model_segments]):
        print("Error: At least one input file must be provided")
        sys.exit(1)
    
    # Check that provided files exist
    for file_arg, file_path in [('copy-ratio', args.copy_ratio), 
                               ('segments', args.segments), 
                               ('model-segments', args.model_segments)]:
        if file_path and not Path(file_path).exists():
            print(f"Error: {file_arg} file not found: {file_path}")
            sys.exit(1)
    
    result = generate_report(args.copy_ratio, args.segments, args.model_segments, 
                           args.sample, args.output)
    
    return result

if __name__ == '__main__':
    main()