#!/usr/bin/env python3
"""
Map gene coordinates to copy number calls from GATK CNV results.

Takes a GTF file with gene coordinates and a GATK calls.cns.tsv file,
and produces a BED file with each gene's copy number call.
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path
import re

def read_gtf_genes(gtf_file):
    """
    Read GTF file and extract gene coordinates.
    Returns DataFrame with columns: chr, start, end, gene_id, gene_name
    """
    genes = []
    
    print(f"Reading GTF file: {gtf_file}")
    
    with open(gtf_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            # Only process gene features
            if parts[2] != 'gene':
                continue
                
            chr_name = parts[0]
            start = int(parts[3]) - 1  # Convert to 0-based coordinates for BED
            end = int(parts[4])
            attributes = parts[8]
            
            # Parse gene_id and gene_name from attributes
            gene_id = None
            gene_name = None
            
            # Look for gene_id
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            if gene_id_match:
                gene_id = gene_id_match.group(1)
            
            # Look for gene_name
            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
            if gene_name_match:
                gene_name = gene_name_match.group(1)
            
            if gene_id:  # Only include if we found a gene_id
                genes.append({
                    'chr': chr_name,
                    'start': start,
                    'end': end,
                    'gene_id': gene_id,
                    'gene_name': gene_name if gene_name else gene_id
                })
            
            if line_num % 100000 == 0:
                print(f"  Processed {line_num} lines, found {len(genes)} genes...")
    
    genes_df = pd.DataFrame(genes)
    print(f"Found {len(genes_df)} genes in GTF file")
    
    return genes_df

def read_gatk_calls(calls_file):
    """
    Read GATK calls.cns.tsv file, skipping SAM headers
    Returns DataFrame with copy number segments
    """
    with open(calls_file, 'r') as f:
        lines = []
        for line in f:
            if not line.startswith('@'):
                lines.append(line)
    
    from io import StringIO
    calls_df = pd.read_csv(StringIO(''.join(lines)), sep='\t')
    
    print(f"Read {len(calls_df)} copy number segments from {calls_file}")
    return calls_df

def map_genes_to_calls(genes_df, calls_df):
    """
    Map each gene to its copy number call based on genomic overlap.
    """
    print("Mapping genes to copy number calls...")
    
    gene_calls = []
    
    for idx, gene in genes_df.iterrows():
        gene_chr = gene['chr']
        gene_start = gene['start']
        gene_end = gene['end']
        
        # Find overlapping segments
        chr_segments = calls_df[calls_df['CONTIG'] == gene_chr]
        
        overlapping_segments = chr_segments[
            (chr_segments['START'] <= gene_end) & 
            (chr_segments['END'] >= gene_start)
        ]
        
        if len(overlapping_segments) == 0:
            # No overlapping segment found
            call = 'N/A'
            mean_log2cr = 'N/A'
        elif len(overlapping_segments) == 1:
            # Single overlapping segment
            seg = overlapping_segments.iloc[0]
            call = seg['CALL']
            mean_log2cr = seg['MEAN_LOG2_COPY_RATIO']
        else:
            # Multiple overlapping segments - choose by maximum overlap
            best_overlap = 0
            best_seg = None
            
            for _, seg in overlapping_segments.iterrows():
                overlap_start = max(gene_start, seg['START'])
                overlap_end = min(gene_end, seg['END'])
                overlap_length = max(0, overlap_end - overlap_start)
                
                if overlap_length > best_overlap:
                    best_overlap = overlap_length
                    best_seg = seg
            
            if best_seg is not None:
                call = best_seg['CALL']
                mean_log2cr = best_seg['MEAN_LOG2_COPY_RATIO']
            else:
                call = 'N/A'
                mean_log2cr = 'N/A'
        
        # Create description
        if call != 'N/A':
            if call == '0':
                call_desc = "neutral"
            elif call == '+':
                call_desc = "amplified"
            elif call == '-':
                call_desc = "deleted"
            else:
                call_desc = f"call_{call}"
            
            description = f"{gene['gene_name']};{call_desc};log2cr={mean_log2cr:.3f}" if mean_log2cr != 'N/A' else f"{gene['gene_name']};{call_desc}"
        else:
            description = f"{gene['gene_name']};no_coverage"
        
        gene_calls.append({
            'chr': gene['chr'],
            'start': gene['start'],
            'end': gene['end'],
            'gene_id': gene['gene_id'],
            'gene_name': gene['gene_name'],
            'call': call,
            'mean_log2cr': mean_log2cr,
            'description': description
        })
        
        if (idx + 1) % 5000 == 0:
            print(f"  Processed {idx + 1}/{len(genes_df)} genes...")
    
    result_df = pd.DataFrame(gene_calls)
    
    # Print summary
    if len(result_df) > 0:
        call_counts = result_df['call'].value_counts()
        print(f"\nCopy number call summary:")
        for call, count in call_counts.items():
            if call == '0':
                print(f"  Neutral: {count} genes")
            elif call == '+':
                print(f"  Amplified: {count} genes")
            elif call == '-':
                print(f"  Deleted: {count} genes")
            elif call == 'N/A':
                print(f"  No coverage: {count} genes")
            else:
                print(f"  {call}: {count} genes")
    
    return result_df

def write_vcf_file(gene_calls_df, output_file, sample_id):
    """
    Write results to VCF format file.
    Uses structural variant format for gene-level CNV calls.
    """
    print(f"Writing VCF file: {output_file}")
    
    from datetime import datetime
    
    with open(output_file, 'w') as f:
        # Write VCF header
        f.write("##fileformat=VCFv4.3\n")
        f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        f.write("##source=GATK_CNV_GeneCalls\n")
        f.write("##reference=GRCh38\n")
        f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
        f.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n")
        f.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
        f.write("##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene name\">\n")
        f.write("##INFO=<ID=GENE_ID,Number=1,Type=String,Description=\"Gene identifier\">\n")
        f.write("##INFO=<ID=LOG2CR,Number=1,Type=Float,Description=\"Log2 copy ratio\">\n")
        f.write("##INFO=<ID=CNV_CALL,Number=1,Type=String,Description=\"Copy number variant call (AMP/DEL/NEU)\">\n")
        f.write("##ALT=<ID=CNV,Description=\"Copy number variant\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n")
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n")
        
        # Only write entries for genes with actual CNV calls (not neutral or no coverage)
        cnv_genes = gene_calls_df[gene_calls_df['call'].isin(['+', '-'])]
        
        for _, gene in cnv_genes.iterrows():
            chrom = gene['chr']
            pos = gene['start'] + 1  # Convert to 1-based coordinates
            gene_id = f"{gene['gene_name']}_{chrom}_{pos}"
            ref = "N"  # Symbolic for structural variants
            alt = "<CNV>"
            qual = "."
            filter_field = "PASS"
            
            # Determine SVTYPE and copy number
            if gene['call'] == '+':
                svtype = "DUP"
                cn = "3"  # Assuming amplification means 3+ copies
                cnv_call = "AMP"
            elif gene['call'] == '-':
                svtype = "DEL" 
                cn = "1"  # Assuming deletion means 1 copy
                cnv_call = "DEL"
            else:
                continue  # Skip neutral/no coverage
            
            # Calculate SVLEN
            svlen = gene['end'] - gene['start']
            
            # Build INFO field
            info_parts = [
                f"SVTYPE={svtype}",
                f"END={gene['end']}",
                f"SVLEN={svlen}",
                f"GENE={gene['gene_name']}",
                f"GENE_ID={gene['gene_id']}",
                f"CNV_CALL={cnv_call}"
            ]
            
            if gene['mean_log2cr'] != 'N/A':
                info_parts.append(f"LOG2CR={gene['mean_log2cr']:.3f}")
            
            info = ";".join(info_parts)
            
            # FORMAT and sample data
            format_field = "GT:CN"
            if gene['call'] == '+':
                sample_data = "./.:3"  # Unknown genotype, 3 copies
            else:  # deletion
                sample_data = "./.:1"  # Unknown genotype, 1 copy
            
            f.write(f"{chrom}\t{pos}\t{gene_id}\tN\t<CNV>\t.\tPASS\t{info}\t{format_field}\t{sample_data}\n")
    
    print(f"Wrote {len(cnv_genes)} CNV calls to {output_file}")
    
    return len(cnv_genes)

def main():
    parser = argparse.ArgumentParser(description='Map gene coordinates to GATK CNV calls')
    parser.add_argument('--gtf', required=True, help='GTF file with gene coordinates')
    parser.add_argument('--calls', required=True, help='GATK calls.cns.tsv file')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--output', required=True, help='Output VCF file')
    
    args = parser.parse_args()
    
    # Validate input files
    if not Path(args.gtf).exists():
        print(f"Error: GTF file not found: {args.gtf}")
        sys.exit(1)
    
    if not Path(args.calls).exists():
        print(f"Error: Calls file not found: {args.calls}")
        sys.exit(1)
    
    print(f"Processing sample: {args.sample}")
    print("=" * 50)
    
    # Read input files
    genes_df = read_gtf_genes(args.gtf)
    calls_df = read_gatk_calls(args.calls)
    
    # Map genes to calls
    gene_calls_df = map_genes_to_calls(genes_df, calls_df)
    
    # Write output
    write_vcf_file(gene_calls_df, args.output, args.sample)
    
    print(f"\nCompleted processing for sample {args.sample}")

if __name__ == '__main__':
    main()