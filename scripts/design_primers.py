#!/usr/bin/env python3
"""
Design primers from multiple sequence alignments using primer3.

This script:
1. Reads a MAFFT alignment file
2. Generates a consensus sequence with conservation scores
3. Identifies conserved regions suitable for primer binding
4. Uses primer3 to design primers in those regions
"""

import argparse
import subprocess
import json
from pathlib import Path
from collections import Counter

# IUPAC ambiguity codes
IUPAC = {
    frozenset(['A']): 'A', frozenset(['C']): 'C', frozenset(['G']): 'G', frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y', frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W', frozenset(['G', 'T']): 'K', frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B', frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H', frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N'
}


def parse_fasta(filepath):
    """Parse FASTA file, return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())
    
    if current_header:
        sequences.append((current_header, ''.join(current_seq)))
    
    return sequences


def generate_consensus(sequences, min_coverage=0.5, conservation_threshold=0.7):
    """
    Generate consensus sequence from aligned sequences.
    Returns: (consensus_sequence, conservation_scores)
    """
    if not sequences:
        return "", []

    # Get alignment length (use max to handle variable length)
    seq_length = max(len(seq) for _, seq in sequences)
    consensus = []
    conservation = []

    for pos in range(seq_length):
        # Count bases at this position (excluding gaps)
        bases = [seq[pos] for _, seq in sequences if pos < len(seq)]
        bases_no_gap = [b for b in bases if b not in '-N']
        
        if not bases_no_gap:
            # All gaps at this position
            consensus.append('-')
            conservation.append(0.0)
            continue
        
        # Calculate frequencies
        counts = Counter(bases_no_gap)
        total = sum(counts.values())
        
        # Get most common base
        most_common_base, most_common_count = counts.most_common(1)[0]
        freq = most_common_count / total
        conservation.append(freq)
        
        # Determine consensus character - always use most common base
        # (primer3 doesn't handle IUPAC well)
        consensus.append(most_common_base)
    
    return ''.join(consensus), conservation


def find_conserved_regions(conservation_scores, min_length=18, threshold=0.8, window=20):
    """Find regions with high conservation scores suitable for primers."""
    regions = []
    n = len(conservation_scores)
    
    i = 0
    while i < n - min_length:
        # Check if this window has high conservation
        window_scores = conservation_scores[i:i+window]
        avg_conservation = sum(window_scores) / len(window_scores)
        
        if avg_conservation >= threshold:
            # Extend the region as long as conservation stays high
            end = i + window
            while end < n:
                next_window = conservation_scores[end:end+10]
                if not next_window or sum(next_window)/len(next_window) < threshold:
                    break
                end += 10
            
            regions.append({
                'start': i,
                'end': end,
                'length': end - i,
                'avg_conservation': sum(conservation_scores[i:end]) / (end - i)
            })
            i = end
        else:
            i += 1
    
    return regions


def run_primer3(sequence, taxon_name, target_regions=None, product_range=(200, 800)):
    """Run primer3 to design primers."""
    # Remove gaps from sequence for primer3
    clean_seq = sequence.replace('-', '').replace('N', '')
    
    if len(clean_seq) < product_range[0]:
        return {'error': f'Sequence too short ({len(clean_seq)} bp)'}
    
    # Build primer3 input
    input_lines = [
        f"SEQUENCE_ID={taxon_name}",
        f"SEQUENCE_TEMPLATE={clean_seq}",
        "PRIMER_TASK=generic",
        "PRIMER_PICK_LEFT_PRIMER=1",
        "PRIMER_PICK_RIGHT_PRIMER=1",
        "PRIMER_OPT_SIZE=20",
        "PRIMER_MIN_SIZE=18",
        "PRIMER_MAX_SIZE=25",
        f"PRIMER_PRODUCT_SIZE_RANGE={product_range[0]}-{product_range[1]}",
        "PRIMER_OPT_TM=60.0",
        "PRIMER_MIN_TM=55.0",
        "PRIMER_MAX_TM=65.0",
        "PRIMER_MIN_GC=40.0",
        "PRIMER_MAX_GC=60.0",
        "PRIMER_NUM_RETURN=5",
        "="
    ]
    
    input_text = '\n'.join(input_lines)
    
    result = subprocess.run(
        ['primer3_core'],
        input=input_text,
        capture_output=True,
        text=True
    )
    
    return parse_primer3_output(result.stdout)


def parse_primer3_output(output):
    """Parse primer3 output into structured format."""
    results = {'primers': [], 'raw': output}

    lines = output.strip().split('\n')
    data = {}
    for line in lines:
        if '=' in line:
            key, value = line.split('=', 1)
            data[key] = value

    num_returned = int(data.get('PRIMER_PAIR_NUM_RETURNED', 0))

    for i in range(num_returned):
        primer = {
            'pair_id': i,
            'left_sequence': data.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''),
            'right_sequence': data.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
            'left_tm': float(data.get(f'PRIMER_LEFT_{i}_TM', 0)),
            'right_tm': float(data.get(f'PRIMER_RIGHT_{i}_TM', 0)),
            'left_gc': float(data.get(f'PRIMER_LEFT_{i}_GC_PERCENT', 0)),
            'right_gc': float(data.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 0)),
            'product_size': int(data.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0)),
            'penalty': float(data.get(f'PRIMER_PAIR_{i}_PENALTY', 999)),
        }
        results['primers'].append(primer)

    return results


def calculate_conservation_stats(conservation_scores):
    """Calculate statistics about conservation across the alignment."""
    if not conservation_scores:
        return {}

    # Filter out gap positions (conservation = 0)
    non_gap = [c for c in conservation_scores if c > 0]
    if not non_gap:
        return {}

    avg_conservation = sum(non_gap) / len(non_gap)
    high_conservation = sum(1 for c in non_gap if c >= 0.9)
    medium_conservation = sum(1 for c in non_gap if 0.7 <= c < 0.9)
    low_conservation = sum(1 for c in non_gap if c < 0.7)

    return {
        'total_positions': len(conservation_scores),
        'non_gap_positions': len(non_gap),
        'avg_conservation': round(avg_conservation, 3),
        'high_conservation_pct': round(100 * high_conservation / len(non_gap), 1),
        'medium_conservation_pct': round(100 * medium_conservation / len(non_gap), 1),
        'low_conservation_pct': round(100 * low_conservation / len(non_gap), 1),
    }


def design_primers_for_taxon(alignment_file, output_dir, taxon_name=None, max_seqs=500):
    """Main function to design primers for a single taxon alignment."""
    alignment_file = Path(alignment_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not taxon_name:
        taxon_name = alignment_file.stem.replace('_aligned', '')

    # Initialize stats tracking
    stats = {
        'taxon': taxon_name,
        'alignment': {},
        'consensus': {},
        'conservation': {},
        'primers': {}
    }

    print(f"\n{'='*60}")
    print(f"🧬 Designing primers for: {taxon_name}")
    print(f"{'='*60}")

    # Step 1: Read alignment
    print("   → Reading alignment...")
    sequences = parse_fasta(alignment_file)
    original_count = len(sequences)
    stats['alignment']['total_sequences'] = original_count
    print(f"     📊 Found {original_count} sequences in alignment")

    if len(sequences) == 0:
        return {'error': 'No sequences in alignment', 'taxon': taxon_name, 'stats': stats}

    # Sample if too many sequences (for speed)
    sampled = False
    if len(sequences) > max_seqs:
        import random
        sequences = random.sample(sequences, max_seqs)
        sampled = True
        stats['alignment']['sampled'] = True
        stats['alignment']['sample_size'] = max_seqs
        print(f"     📊 Sampled {max_seqs} sequences for consensus (from {original_count})")
    else:
        stats['alignment']['sampled'] = False

    # Get alignment length
    if sequences:
        align_length = len(sequences[0][1])
        stats['alignment']['alignment_length'] = align_length
        print(f"     📊 Alignment length: {align_length} bp")

    # Step 2: Generate consensus
    print("   → Generating consensus sequence...")
    consensus, conservation = generate_consensus(sequences)
    clean_consensus = consensus.replace('-', '')

    stats['consensus']['raw_length'] = len(consensus)
    stats['consensus']['clean_length'] = len(clean_consensus)
    stats['consensus']['gap_positions'] = len(consensus) - len(clean_consensus)

    print(f"     📊 Consensus length: {len(clean_consensus)} bp (gaps removed: {stats['consensus']['gap_positions']})")

    # Save consensus
    consensus_file = output_dir / f"{taxon_name}_consensus.fasta"
    with open(consensus_file, 'w') as f:
        f.write(f">{taxon_name}_consensus\n{clean_consensus}\n")

    # Step 3: Find conserved regions and calculate conservation stats
    print("   → Analyzing conservation...")
    conservation_stats = calculate_conservation_stats(conservation)
    stats['conservation'] = conservation_stats

    if conservation_stats:
        print(f"     📊 Average conservation: {conservation_stats['avg_conservation']:.1%}")
        print(f"        High (≥90%): {conservation_stats['high_conservation_pct']:.1f}% of positions")
        print(f"        Medium (70-90%): {conservation_stats['medium_conservation_pct']:.1f}% of positions")
        print(f"        Low (<70%): {conservation_stats['low_conservation_pct']:.1f}% of positions")

    print("   → Finding conserved regions...")
    conserved = find_conserved_regions(conservation)
    stats['conservation']['num_conserved_regions'] = len(conserved)
    print(f"     📊 Found {len(conserved)} conserved regions suitable for primers")

    # Step 4: Design primers
    print("   → Running primer3...")
    primer_results = run_primer3(consensus, taxon_name)

    if 'error' in primer_results:
        print(f"     ⚠ {primer_results['error']}")
        stats['primers']['error'] = primer_results['error']
        return {'error': primer_results['error'], 'taxon': taxon_name, 'stats': stats}

    stats['primers']['num_pairs'] = len(primer_results.get('primers', []))

    if primer_results['primers']:
        print(f"     ✓ Designed {len(primer_results['primers'])} primer pairs")

        # Add stats to primer results
        primer_results['stats'] = stats

        # Save results
        results_file = output_dir / f"{taxon_name}_primers.json"
        with open(results_file, 'w') as f:
            json.dump(primer_results, f, indent=2)

        # Print best primer pair
        best = primer_results['primers'][0]
        stats['primers']['best_forward'] = best['left_sequence']
        stats['primers']['best_reverse'] = best['right_sequence']
        stats['primers']['best_product_size'] = best['product_size']

        print(f"\n   📌 Best primer pair:")
        print(f"      Forward: {best['left_sequence']} (Tm={best['left_tm']:.1f}°C, GC={best['left_gc']:.0f}%)")
        print(f"      Reverse: {best['right_sequence']} (Tm={best['right_tm']:.1f}°C, GC={best['right_gc']:.0f}%)")
        print(f"      Product: {best['product_size']} bp")
    else:
        print("     ⚠ No suitable primers found")
        # Save stats even if no primers found
        stats_file = output_dir / f"{taxon_name}_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)

    return {
        'taxon': taxon_name,
        'num_sequences': original_count,
        'sequences_used': len(sequences),
        'consensus_length': len(clean_consensus),
        'conserved_regions': len(conserved),
        'primers': primer_results.get('primers', []),
        'stats': stats
    }


def main():
    parser = argparse.ArgumentParser(description='Design primers from alignment')
    parser.add_argument('--alignment', '-a', help='Path to alignment FASTA file')
    parser.add_argument('--taxon', '-t', help='Taxon name (defaults to filename)')
    parser.add_argument('--output', '-o', help='Output directory', default='.')
    parser.add_argument('--gene', '-g', choices=['18S', '16S'], default='18S')
    parser.add_argument('--all', action='store_true', help='Process all alignments')
    parser.add_argument('--max-seqs', type=int, default=500, help='Max sequences for consensus')
    args = parser.parse_args()

    if args.all:
        # Process all alignments in the gene directory
        base_dir = Path(__file__).parent.parent / args.gene
        alignments = list(base_dir.glob('*/align/*_aligned.fasta'))

        print(f"🧬 Primer Design Pipeline - {args.gene}")
        print(f"Found {len(alignments)} alignments to process\n")

        results = []
        for align_file in alignments:
            if align_file.stat().st_size == 0:
                print(f"Skipping empty: {align_file.parent.parent.name}")
                continue

            output_dir = align_file.parent.parent / 'primers'
            result = design_primers_for_taxon(align_file, output_dir, max_seqs=args.max_seqs)
            results.append(result)

        # Summary
        print(f"\n{'='*60}")
        print("📊 Summary")
        print(f"{'='*60}")
        success = [r for r in results if r.get('primers')]
        print(f"Successfully designed primers for {len(success)}/{len(results)} taxa")

    elif args.alignment:
        output_dir = args.output or Path(args.alignment).parent.parent / 'primers'
        design_primers_for_taxon(args.alignment, output_dir, args.taxon, args.max_seqs)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()

