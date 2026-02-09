# Methods: Paramoeba Pipeline

## Overview

This document describes how the Paramoeba sequences were processed through the primer design pipeline.

## Step 1: Sequence Extraction

**Command:**
```bash
python scripts/extract_sequences_by_taxon.py \
    --tsv eukcensus_metadata/eukcensus_18S.clusters.97.tsv \
    --taxon Paramoeba \
    --rank genus \
    --output-dir primer_design_pipeline/18S/Paramoeba \
    --fasta eukcensus_metadata/eukcensus_2025_18S.97.fna \
    --taxon-original Paramoeba_G \
    --gene 18S
```

**Result:** 20 sequences extracted

## Step 2: Length Filtering (Optional)

**Command:**
```bash
seqkit seq -m 1200 Paramoeba_18S.fasta > filtered/Paramoeba_1200bp.fasta
```

**Note:** This step was skipped because seqkit was not installed. All 20 original sequences were used.

## Step 3: Multiple Sequence Alignment

**Command:**
```bash
mafft --auto --thread 4 Paramoeba_18S.fasta > align/Paramoeba_aligned.fasta
```

**Settings:**
- Algorithm: Auto-selected by MAFFT based on sequence count
- Threads: 4

**Result:** 20 sequences aligned

## Step 4: Primer Design

See `primers/PRIMER_METHODS.md` for primer design details.

## Output Files

| File | Description |
|------|-------------|
| `Paramoeba_18S.fasta` | Original extracted sequences |
| `Paramoeba_sequences.tsv` | Sequence IDs and metadata |
| `align/Paramoeba_aligned.fasta` | MAFFT alignment |
| `primers/Paramoeba_consensus.fasta` | Consensus sequence |
| `primers/Paramoeba_primers.json` | Designed primers |

