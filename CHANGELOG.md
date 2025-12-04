# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] - 2025-12-04

### Added
- Complete lysozyme annotation pipeline with 7-step workflow
- BLAST-based homology search (tblastn)
- Quality filtering with identity and score thresholds
- Optional SSEARCH Smith-Waterman realignment
- BEDTools region merging
- Score density calculation for best-hit selection
- Comprehensive pseudogene detection with 6 mutation types
- Batch processing for multiple genomes
- Comparative analysis with statistical tests
- GFF3 export functionality
- Automatic dependency installation (BLAST+, BEDTools, FASTA36)
- Multi-threading support for parallelization
- Detailed logging with adjustable verbosity
- Professional documentation (README, examples)

### Features
- **Mutation Detection:**
  - Non-synonymous substitutions
  - In-frame indels
  - Frameshifts
  - Premature stop codons
  - Missing start/stop codons (genomic verification)

- **Comparative Analysis:**
  - Pathogenic vs non-pathogenic comparison
  - Protein presence/absence matrix
  - Core/accessory/unique protein classification
  - Differential protein identification
  - Statistical tests (Mann-Whitney U)
  - Hierarchical clustering and PCA
  - Publication-quality visualizations

### Optimizations
- Smart start/stop codon verification (alignment check before genomic search)
- Minimal logging verbosity (3-5 lines per genome)
- Efficient multi-threading for BLAST and SSEARCH

### Technical Details
- Python 3.8+
- BioPython, pandas, scipy, matplotlib, seaborn, scikit-learn
- BLAST+ 2.10.0+, BEDTools 2.29.0+, FASTA36 (optional)
