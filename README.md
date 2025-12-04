# Lysozyme Annotation Pipeline

Automated pipeline for homology search, filtering, realignment, and pseudogene characterization of lysozyme genes in bacterial genomes.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Single Genome Analysis](#single-genome-analysis)
  - [Batch Processing](#batch-processing)
  - [Comparative Analysis](#comparative-analysis)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Methods](#methods)
- [Parameters](#parameters)
- [Troubleshooting](#troubleshooting)

## Overview

This pipeline implements a workflow for identifying and characterizing lysozyme genes and pseudogenes in bacterial genomes through homology analysis and inactivating mutation detection.

### Key Features

- BLAST-based homology search (tblastn) with optimized parameters
- Quality filtering using strict identity and score thresholds
- Rigorous realignment with SSEARCH (Smith-Waterman algorithm) - **Mandatory**
- Adjacent HSP merging using BEDTools
- Best-hit selection via score density calculation
- Comprehensive pseudogene detection and mutation quantification
- Batch processing for multiple genomes
- Advanced comparative analysis for pathogenic vs non-pathogenic comparison

## Requirements

### External Tools

- NCBI BLAST+ (>= 2.10.0) - **Required**
- BEDTools (>= 2.29.0) - **Required**
- FASTA36 suite (ssearch36) - **Required**

The pipeline can automatically download and install these tools locally if not found in system PATH.

### Python Dependencies

- Python >= 3.8
- BioPython >= 1.79
- pandas >= 1.3.0
- scipy >= 1.7.0 (for comparative analysis)
- matplotlib >= 3.4.0 (for visualizations)
- seaborn >= 0.11.0 (for visualizations)
- scikit-learn >= 0.24.0 (for PCA)

## Installation

1. Clone or download the project:
```bash
cd ~/Documents/lysozyme_pipeline
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

3. Verify external tools (optional - pipeline will install if missing):
```bash
makeblastdb -version
tblastn -version
bedtools --version
ssearch36 -h
```

## Usage

### Single Genome Analysis

Basic usage for analyzing a single genome:

```bash
python pipeline.py \
    --input genome.fasta \
    --lysozymes reference_lysozymes.fasta \
    --output results/ \
    --num-threads 8
```

### Batch Processing

Process multiple genomes in parallel:

```bash
python pipeline.py \
    --input-dir genomes/ \
    --lysozymes reference_lysozymes.fasta \
    --output batch_results/ \
    --num-threads 20
```

This will:
- Process each genome independently
- Generate per-genome results in subdirectories
- Create aggregated comparison table (`all_genomes_aggregated.tsv`)
- Generate summary statistics across all genomes

### Comparative Analysis

After batch processing, perform statistical comparison between pathogenic and non-pathogenic groups:

1. Create metadata file (`metadata.tsv`):
```tsv
genome_id	pathogenicity	species	strain
genome1	pathogenic	E. coli	O157:H7
genome2	non-pathogenic	E. coli	K-12
genome3	pathogenic	Salmonella enterica	Typhimurium
```

2. Run comparative analysis:
```bash
python comparative_analysis.py \
    --batch-output batch_results/ \
    --metadata metadata.tsv \
    --output-dir comparative_results/
```

This generates:
- Statistical tests (Mann-Whitney U) for pseudogene rates, score density, frameshifts
- Protein classification (core/accessory/unique)
- Differential protein identification
- Five publication-quality visualizations (heatmaps, PCA, dendrogram, boxplots)
- Executive report (Markdown) with automated conclusions

## Pipeline Steps

The pipeline implements a seven-step workflow:

### 1. BLAST Database Preparation
Creates nucleotide database from genome FASTA using `makeblastdb`.

### 2. tblastn Search
Searches protein queries against nucleotide database translated in all six reading frames.

Parameters:
- E-value threshold: 1e-10
- Matrix: BLOSUM62
- Output format: tabular with sequences

### 3. Quality Filtering
Removes low-quality alignments based on:
- Minimum identity: 20% (default)
- Minimum BLOSUM62 score: 120 (default)
- Discard rule: BOTH conditions must fail to remove

### 4. SSEARCH Realignment
Re-aligns filtered hits using Smith-Waterman algorithm with statistical permutation testing.

Parameters:
- Shuffles: 1000 permutations
- Statistical significance threshold: p < 0.001
- **Mandatory step** for accurate pseudogene detection

### 5. Region Merging
Merges overlapping/adjacent HSPs using BEDTools merge.

Parameters:
- Maximum distance: 50 bp
- Strand-aware merging

### 6. Score Density Calculation
For each merged region:
- Groups HSPs by protein ID
- Calculates score density = total_score / total_length
- Selects best protein (highest score density)

### 7. Pseudogene Detection
Identifies inactivating mutations:
- Non-synonymous substitutions (query != subject amino acid)
- In-frame indels (gaps multiple of 3)
- Frameshifts (gaps not multiple of 3)
- Premature stop codons (internal * in alignment)
- Missing start codon (verified in genomic DNA)
- Missing stop codon (verified in genomic DNA)

Start/stop codon verification strategy:
1. Check if alignment includes start (M at position 1) or stop (reaches query end)
2. If not, search ±150bp in genomic DNA
3. Mark as missing if not found

## Output Files

### Single Genome Mode

```
output/
├── blast/
│   ├── genome_db.*                    # BLAST database files
│   └── blast_results.tsv              # Raw BLAST output
├── ssearch/
│   └── ssearch_results.tsv            # SSEARCH realignments
├── merge/
│   ├── blast_hits.bed                 # BED format regions
│   ├── merged_regions.bed             # Merged regions
│   └── merged_regions.tsv             # Merged regions table
└── final/
    ├── region_annotations.tsv         # Score density annotations
    ├── pseudogene_annotations.tsv     # Mutation analysis
    └── lysozyme_annotations.gff3      # GFF3 format output
```

### Batch Mode

```
batch_results/
├── genome1/
│   ├── blast/...
│   ├── merge/...
│   └── final/...
├── genome2/...
├── comparative_analysis/
│   ├── all_genomes_aggregated.tsv     # Combined results
│   ├── summary_stats.csv              # Per-genome summary
│   ├── presence_absence_matrix.tsv    # Protein × genome matrix
│   └── comparative_summary.txt        # Text report
└── batch_gff3/
    ├── genome1.gff3
    └── genome2.gff3
```

### Comparative Analysis Mode

```
comparative_results/
├── EXECUTIVE_REPORT.md                # Automated report with conclusions
├── protein_profiles.tsv               # Per-protein statistics
├── protein_categories.txt             # Core/accessory/unique classification
├── group_statistics.tsv               # Pathogenic vs non-pathogenic stats
├── differential_proteins.tsv          # Enriched proteins
├── pathogenicity_differential.txt     # Detailed report
├── presence_absence_matrix.tsv        # Binary matrix
├── pseudogenization_matrix.tsv        # Pseudogenization rate matrix
└── figures/
    ├── heatmap_presence_absence.png   # Protein × genome heatmap
    ├── heatmap_pseudogenization.png   # Pseudogenization heatmap
    ├── pca_lysozyme_profiles.png      # PCA with variance explained
    ├── dendrogram_similarity.png      # Hierarchical clustering
    └── boxplots_group_comparison.png  # Statistical comparisons
```

## Methods

### Score Density Algorithm

For each merged genomic region containing multiple protein hits:

```
For each protein P with HSPs {h1, h2, ..., hn}:
    total_score = sum(hsp.score for hsp in HSPs)
    total_length = sum(hsp.length for hsp in HSPs)
    score_density = total_score / total_length

best_protein = protein with highest score_density
```

This metric balances alignment quality (score) with coverage (length), providing more reliable protein identification than raw score alone.

### Pseudogene Classification

A region is classified as a pseudogene if:

```
total_disablements = non_synonymous + in_frame_indels + frameshifts + 
                     premature_stops + missing_start + missing_stop

if total_disablements >= min_disablements (default: 1):
    is_pseudogene = True
```

### Statistical Tests (Comparative Analysis)

Mann-Whitney U test (non-parametric) for:
- Pseudogene rate per genome
- Mean score density
- Total frameshifts per genome

Differential protein identification:
- Presence difference >= 20% between groups
- OR pseudogenization rate difference >= 20%

Core protein definition:
- Present in >= 90% of genomes

## Parameters

### BLAST Parameters
- `--evalue`: E-value threshold (default: 1e-10)
- `--num-threads`: CPU threads for parallelization (default: 1)

### Filtering Parameters
- `--min-identity`: Minimum % identity (default: 20)
- `--min-score`: Minimum BLOSUM62 score (default: 120)

### SSEARCH Parameters
- `--ssearch-shuffles`: Number of permutations (default: 1000)

### Pseudogene Detection
- `--min-disablements`: Minimum mutations for pseudogene classification (default: 1)

### Comparative Analysis
- `--alpha`: Statistical significance level (default: 0.05)

## Troubleshooting

### BLAST+ not found
The pipeline will attempt automatic installation to `BLAST/` directory. If this fails:
```bash
sudo apt install ncbi-blast+  # Debian/Ubuntu
brew install blast            # macOS
```

### BEDTools compilation failed
```bash
sudo apt install build-essential zlib1g-dev
```

### SSEARCH not available
SSEARCH is required for the pipeline. To install:
```bash
sudo apt install fasta3       # Debian/Ubuntu
```
Alternatively, the pipeline will attempt automatic installation.

### Low pseudogene count
Check alignment quality with `--verbose` flag. Consider adjusting:
- `--min-identity` (lower for more divergent sequences)
- `--min-score` (lower threshold accepts more hits)

### High pseudogene count
Verify reference proteins are appropriate for target organism. Consider:
- Using more closely related reference sequences
- Increasing `--min-disablements` threshold

### Comparative analysis errors
Ensure:
- Metadata `genome_id` matches genome filenames exactly
- `pathogenicity` values are: `pathogenic`, `non-pathogenic`, or `unknown`
- TSV file uses tab separators (not spaces)

