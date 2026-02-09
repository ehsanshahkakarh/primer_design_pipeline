# Methods: Vannellidae Pipeline

## Overview

This document describes how the Vannellidae sequences were processed through the primer design pipeline.

## Step 1: Sequence Extraction

**Command:**
```bash
python scripts/extract_sequences_by_taxon.py \
    --tsv eukcensus_metadata/eukcensus_18S.clusters.97.tsv \
    --taxon Vannellidae \
    --rank family \
    --output-dir primer_design_pipeline/18S/Vannellidae \
    --fasta eukcensus_metadata/eukcensus_2025_18S.97.fna \
    --taxon-original Vannellidae_F \
    --gene 18S
```

**Result:** 56 sequences extracted

## Step 2: Length Filtering (Optional)

**Command:**
```bash
seqkit seq -m 1200 Vannellidae_18S.fasta > filtered/Vannellidae_1200bp.fasta
```

**Note:** This step was skipped because seqkit was not installed. All 56 original sequences were used.

## Step 3: Multiple Sequence Alignment

**Command:**
```bash
mafft --auto --thread 4 Vannellidae_18S.fasta > align/Vannellidae_aligned.fasta
```

**Settings:**
- Algorithm: Auto-selected by MAFFT based on sequence count
- Threads: 4

**Result:** 56 sequences aligned

## Step 4: Primer Design

See `primers/PRIMER_METHODS.md` for primer design details.

## Output Files

| File | Description |
|------|-------------|
| `Vannellidae_18S.fasta` | Original extracted sequences |
| `Vannellidae_sequences.tsv` | Sequence IDs and metadata |
| `align/Vannellidae_aligned.fasta` | MAFFT alignment |
| `primers/Vannellidae_consensus.fasta` | Consensus sequence |
| `primers/Vannellidae_primers.json` | Designed primers |

