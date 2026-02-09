# Primer Design Methods: Vannella

## Overview

Primers were designed from the multiple sequence alignment using consensus-based primer design with primer3.

## Step 1: Consensus Generation

**Method:** Majority-rule consensus
- At each position, the most common nucleotide was selected
- Gap positions were removed from the final consensus
- No IUPAC ambiguity codes used (primer3 compatibility)

**Input:** 56 aligned sequences
**Output:** Consensus sequence (see `Vannella_consensus.fasta`)

## Step 2: Primer3 Design

**Command (internal):**
```bash
primer3_core << EOF
SEQUENCE_ID=Vannella
SEQUENCE_TEMPLATE=<consensus_sequence>
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_PRODUCT_SIZE_RANGE=200-800
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=55.0
PRIMER_MAX_TM=65.0
PRIMER_MIN_GC=40.0
PRIMER_MAX_GC=60.0
PRIMER_NUM_RETURN=5
=
EOF
```

## Primer3 Settings

| Parameter | Value | Description |
|-----------|-------|-------------|
| PRIMER_OPT_SIZE | 20 bp | Optimal primer length |
| PRIMER_MIN_SIZE | 18 bp | Minimum primer length |
| PRIMER_MAX_SIZE | 25 bp | Maximum primer length |
| PRIMER_OPT_TM | 60.0°C | Optimal melting temperature |
| PRIMER_MIN_TM | 55.0°C | Minimum melting temperature |
| PRIMER_MAX_TM | 65.0°C | Maximum melting temperature |
| PRIMER_MIN_GC | 40% | Minimum GC content |
| PRIMER_MAX_GC | 60% | Maximum GC content |
| PRIMER_PRODUCT_SIZE_RANGE | 200-800 bp | Target amplicon size |
| PRIMER_NUM_RETURN | 5 | Number of primer pairs to return |

## Results

**Best Primer Pair:**
- Forward: `TGCCGACCAGGGATTAGAGA`
- Reverse: `TCCACCAACTAAGAACGGCC`
- Product Size: 329 bp
- Forward Tm: 60.0°C, GC: 55%
- Reverse Tm: 60.0°C, GC: 55%

**Total Primer Pairs Designed:** 5

See `Vannella_primers.json` for all primer pairs with detailed statistics.

