# Methods: Spizellomycetaceae Pipeline

## Overview

This document describes how the Spizellomycetaceae sequences were processed through the primer design pipeline.

## Step 1: Sequence Extraction

**Command:**
```bash
python scripts/extract_sequences_by_taxon.py \
    --tsv eukcensus_metadata/eukcensus_18S.clusters.97.tsv \
    --taxon Spizellomycetaceae \
    --rank family \
    --output-dir primer_design_pipeline/18S/Spizellomycetaceae \
    --fasta eukcensus_metadata/eukcensus_2025_18S.97.fna \
    --taxon-original Spizellomycetaceae_F \
    --gene 18S
```

**Result:** 52 sequences extracted

## Step 2: Length Filtering (Optional)

**Command:**
```bash
seqkit seq -m 1200 Spizellomycetaceae_18S.fasta > filtered/Spizellomycetaceae_1200bp.fasta
```

**Note:** This step was skipped because seqkit was not installed. All 52 original sequences were used.

## Step 3: Multiple Sequence Alignment

**Command:**
```bash
mafft --auto --thread 4 Spizellomycetaceae_18S.fasta > align/Spizellomycetaceae_aligned.fasta
```

**Settings:**
- Algorithm: Auto-selected by MAFFT based on sequence count
- Threads: 4

**Result:** 52 sequences aligned

## Step 4: Primer Design

See `primers/PRIMER_METHODS.md` for primer design details.

## Output Files

| File | Description |
|------|-------------|
| `Spizellomycetaceae_18S.fasta` | Original extracted sequences |
| `Spizellomycetaceae_sequences.tsv` | Sequence IDs and metadata |
| `align/Spizellomycetaceae_aligned.fasta` | MAFFT alignment |
| `primers/Spizellomycetaceae_consensus.fasta` | Consensus sequence |
| `primers/Spizellomycetaceae_primers.json` | Designed primers |

