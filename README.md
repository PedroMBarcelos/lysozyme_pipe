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

After batch processing, perform statistical comparison between different pathotype groups:

1. Create metadata file (`metadata.tsv`) with pathovar information:
```tsv
genome_id	pathovar	species	strain
562_13551	EHEC	Escherichia coli	O157:H7
83334_298	NOT	Escherichia coli	K-12 MG1655
199310_4	AIEC	Escherichia coli	LF82
574521_7	EHEC	Escherichia coli	O104:H4
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
- Differential protein identification between groups
- Five publication-quality visualizations grouped by pathovar:
  - Presence/absence heatmap
  - Pseudogenization heatmap  
  - PCA with pathovar-based coloring
  - Hierarchical dendrogram with pathovar-colored labels
  - Boxplots comparing metrics across pathovars
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
2. If not, search in genomic DNA with **in-frame search** (step=3) within dynamic padding window
3. Dynamic padding = min(150 bp, 20% of reference protein length × 3)
4. Mark as missing if not found within the padded region

### Coverage Metrics

Three coverage metrics are calculated:

1. **reference_coverage_nt**: Alignment length × 3 (estimated nucleotide coverage)
2. **alignment_coverage_nt**: Exact genomic coverage from subject coordinates (sstart to send)
3. **alignment_genomic_ratio**: alignment_coverage_nt / reference_coverage_nt (compactness measure, ~1.0 for no gaps/frameshifts)

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

### Input/Output Parameters

**Required:**
- `-g, --genome FASTA`: Single genome FASTA file (mutually exclusive with `--input-dir`)
- `--input-dir DIR`: Directory with multiple FASTA files for batch processing (mutually exclusive with `-g`)
- `-l, --lysozymes FASTA`: Reference lysozyme protein sequences (UniProtKB/Swiss-Prot format)
- `-o, --output DIR`: Output directory for all results

### Filtering Parameters

- `--min-identity FLOAT`: Minimum percent identity for BLAST filtering (default: 20.0)
- `--min-score INT`: Minimum BLOSUM62 score for BLAST filtering (default: 120)
- `--min-disablements INT`: Minimum number of inactivating mutations to classify as pseudogene (default: 1)
- `--min-coverage FLOAT`: Minimum alignment coverage ratio for reporting (default: 0.8)
- `--final-min-identity FLOAT`: Final identity filter applied after all processing (default: 0.0 - disabled)

### Execution Options

- `--num-threads INT`: Number of CPU threads for SSEARCH parallel execution (default: CPU_COUNT - 1)
- `-v, --verbose`: Enable detailed debug messages
- `--log-file FILE`: Save logs to specified file (optional)

### Comparative Analysis Parameters

When using `comparative_analysis.py`:
- `--batch-output DIR`: Directory containing batch processing results
- `--metadata FILE`: TSV file with genome metadata (genome_id, pathovar, species, strain)
- `--output-dir DIR`: Output directory for comparative results
- `--alpha FLOAT`: Statistical significance level for tests (default: 0.05)

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
- Metadata `genome_id` matches genome filenames exactly (or genome IDs from BVBRC metadata)
- `pathovar` column exists with values like: `EHEC`, `AIEC`, `NOT`, etc.
- TSV file uses tab separators (not spaces)
- genome_id format consistency: dots (562.13551) are automatically normalized to underscores (562_13551)

### Metadata File Format

The metadata file should be tab-separated (TSV) with these columns:

```tsv
genome_id	pathovar	species	strain
562_13551	EHEC	Escherichia coli	O157:H7
83334_298	NOT	Escherichia coli	K-12 MG1655
199310_4	AIEC	Escherichia coli	LF82
```

**Important:**
- `genome_id`: Must match FASTA filenames (without extension) or BVBRC genome IDs
- `pathovar`: Pathotype/pathovar designation (e.g., EHEC, EPEC, AIEC, NOT for non-pathogenic)
- Visualizations will group genomes by `pathovar` values directly
- Use consistent naming: genome IDs with dots are normalized to underscores automatically

## Genome Download Utility

The pipeline includes `download_genomes.py` to fetch genomes from NCBI using BVBRC metadata CSV files:

```bash
python download_genomes.py BVBRC_genome_combined.csv downloaded_genomes/
```

This will:
- Extract genome IDs and assembly accessions from the CSV
- Download FASTA files from NCBI using efetch
- Create a metadata.tsv file with genome_id and pathovar columns
- Normalize genome IDs (dots → underscores) for pipeline compatibility
- Include all non-pathogenic genomes (marked with `NOT` pathovar)

The combined CSV (`BVBRC_genome_combined.csv`) includes:
- 148 pathogenic/clinical E. coli genomes (EHEC, STEC, UPEC, AIEC, ETEC, etc.)
- 14 non-pathogenic reference strains including:
  - Laboratory strains: K-12 MG1655, BW25113, K-12 J53
  - Commensal isolates: HS, HS4, HS1496, HS13-1, HS30-1, 2 HS-C
  - Clinical non-pathogenic isolates: MSHS 472, MSHS 133, 6535NalR