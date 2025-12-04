#!/usr/bin/env python3
"""
Advanced Comparative Analysis Script
Comparative analysis between pathogenic and non-pathogenic genomes.

Usage:
    python comparative_analysis.py --batch-output output/batch_test \
                                   --metadata metadata.tsv \
                                   --output-dir output/batch_test/advanced_comparative
"""

import argparse
import logging
import sys
from pathlib import Path
import pandas as pd

from src.protein_profile_analysis import (
    calculate_protein_profiles,
    classify_protein_categories,
    get_protein_genome_matrix,
    get_pseudogenization_matrix,
    save_protein_profiles
)
from src.pathogenicity_comparison import (
    validate_metadata,
    check_group_balance,
    calculate_group_statistics,
    perform_statistical_tests,
    identify_differential_proteins,
    save_pathogenicity_analysis
)
from src.visualization import generate_all_visualizations


# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def load_batch_results(batch_output_dir: Path) -> pd.DataFrame:
    """
    Load aggregated batch results.
    
    Args:
        batch_output_dir: Batch output directory
    
    Returns:
        DataFrame with aggregated results
    """
    comparative_dir = batch_output_dir / "comparative_analysis"
    aggregated_file = comparative_dir / "all_genomes_aggregated.tsv"
    
    if not aggregated_file.exists():
        raise FileNotFoundError(
            f"Aggregated file not found: {aggregated_file}\n"
            f"Run the batch pipeline first."
        )
    
    logger.debug(f"Loading results from: {aggregated_file}")
    df = pd.read_csv(aggregated_file, sep='\t')
    logger.info(f"{len(df)} regions loaded from {df['genome_id'].nunique()} genomes")
    
    return df


def load_metadata(metadata_path: Path) -> pd.DataFrame:
    """
    Load metadata file.
    
    Args:
        metadata_path: Path to metadata TSV file
    
    Returns:
        DataFrame with metadata
    """
    logger.debug(f"Loading metadata from: {metadata_path}")
    
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    
    # Ignore comment lines (starting with #)
    metadata_df = pd.read_csv(metadata_path, sep='\t', comment='#')
    logger.info(f"Metadata loaded for {len(metadata_df)} genomes")
    
    return metadata_df


def generate_executive_report(
    aggregated_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    protein_profiles: pd.DataFrame,
    protein_categories: dict,
    group_stats: pd.DataFrame,
    statistical_tests: dict,
    differential_proteins: pd.DataFrame,
    output_dir: Path
) -> None:
    """
    Generate consolidated executive report in Markdown format.
    
    Args:
        aggregated_df: Aggregated DataFrame
        metadata_df: DataFrame with metadata
        protein_profiles: Protein profiles
        protein_categories: Protein categories
        group_stats: Statistics per group
        statistical_tests: Test results
        differential_proteins: Differential proteins
        output_dir: Output directory
    """
    logger.debug("Generating executive report...")
    
    report_path = output_dir / "EXECUTIVE_REPORT.md"
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("# Executive Report - Lysozyme Comparative Analysis\n\n")
        f.write("---\n\n")
        
        # SECTION 1: Overview
        f.write("## 1. Overview\n\n")
        
        total_genomes = aggregated_df['genome_id'].nunique()
        n_pathogenic = len(metadata_df[metadata_df['pathogenicity'] == 'pathogenic'])
        n_non_pathogenic = len(metadata_df[metadata_df['pathogenicity'] == 'non-pathogenic'])
        total_regions = len(aggregated_df)
        total_proteins = aggregated_df['best_protein_id'].nunique()
        
        f.write(f"- **Total genomes analyzed:** {total_genomes}\n")
        f.write(f"  - Pathogenic: {n_pathogenic}\n")
        f.write(f"  - Non-pathogenic: {n_non_pathogenic}\n")
        f.write(f"- **Total lysozyme regions identified:** {total_regions}\n")
        f.write(f"- **Unique reference proteins:** {total_proteins}\n\n")
        
        # SECTION 2: Core Lysozyme Repertoire
        f.write("## 2. Core Lysozyme Repertoire\n\n")
        
        f.write(f"**Core lysozymes** (present in ≥90% of genomes):\n\n")
        if len(protein_categories['core']) > 0:
            for protein in protein_categories['core']:
                protein_data = protein_profiles[protein_profiles['protein_id'] == protein].iloc[0]
                f.write(f"- `{protein}`\n")
                f.write(f"  - Presence: {protein_data['presence_rate']:.1f}% of genomes\n")
                f.write(f"  - Pseudogenization rate: {protein_data['pseudogene_rate']:.1f}%\n")
                f.write(f"  - Mean score density: {protein_data['mean_score_density']:.2f}\n\n")
        else:
            f.write("*No core lysozymes identified (threshold=90%)*\n\n")
        
        f.write(f"**Accessory lysozymes:** {len(protein_categories['accessory'])}\n\n")
        f.write(f"**Unique lysozymes:** {len(protein_categories['unique'])}\n\n")
        
        # SECTION 3: Statistical Differences
        f.write("## 3. Pathogenic vs Non-Pathogenic Comparison\n\n")
        
        f.write("### Statistics by Group\n\n")
        f.write("| Metric | Pathogenic | Non-Pathogenic |\n")
        f.write("|---------|-------------|----------------|\n")
        
        path_stats = group_stats[group_stats['group'] == 'pathogenic'].iloc[0] if len(group_stats[group_stats['group'] == 'pathogenic']) > 0 else None
        non_path_stats = group_stats[group_stats['group'] == 'non-pathogenic'].iloc[0] if len(group_stats[group_stats['group'] == 'non-pathogenic']) > 0 else None
        
        if path_stats is not None and non_path_stats is not None:
            f.write(f"| Pseudogene rate | {path_stats['mean_pseudogene_rate']:.2f}% | {non_path_stats['mean_pseudogene_rate']:.2f}% |\n")
            f.write(f"| Mean score density | {path_stats['mean_score_density']:.2f} | {non_path_stats['mean_score_density']:.2f} |\n")
            f.write(f"| Total frameshifts | {path_stats['total_frameshifts']} | {non_path_stats['total_frameshifts']} |\n")
            f.write(f"| Premature stops | {path_stats['total_premature_stops']} | {non_path_stats['total_premature_stops']} |\n")
            f.write(f"| Missing starts | {path_stats['total_missing_start']} | {non_path_stats['total_missing_start']} |\n")
            f.write(f"| Missing stops | {path_stats['total_missing_stop']} | {non_path_stats['total_missing_stop']} |\n\n")
        
        f.write("### Statistical Tests\n\n")
        for test_name, result in statistical_tests.items():
            f.write(f"**{test_name.replace('_', ' ').title()}**\n\n")
            f.write(f"- Teste: {result['test']}\n")
            f.write(f"- p-value: {result['p_value']:.4e}\n")
            f.write(f"- **Result: {'SIGNIFICANT \u2713' if result['significant'] else 'Not significant'}**\n\n")
        
        # SECTION 4: Specific Lysozymes
        f.write("## 4. Differentially Present Lysozymes\n\n")
        
        if not differential_proteins.empty:
            f.write(f"Identified **{len(differential_proteins)} lysozymes** with significant presence or pseudogenization differences between groups.\n\n")
            
            f.write("### Top 10 Differential Lysozymes\n\n")
            f.write("| Protein | Pathogenic Presence | Non-Path Presence | Enriched in |\n")
            f.write("|----------|---------------------|-------------------|-------------|\n")
            
            for _, row in differential_proteins.head(10).iterrows():
                f.write(f"| `{row['protein_id']}` | {row['pathogenic_presence']:.1f}% | {row['non_pathogenic_presence']:.1f}% | {row['enriched_in']} |\n")
            
            f.write("\n")
        else:
            f.write("*No differential lysozymes identified.*\n\n")
        
        # SECTION 5: Visualizations
        f.write("## 5. Visualizations\n\n")
        f.write("Figures generated in `figures/`:\n\n")
        f.write("- `heatmap_presence_absence.png` - Heatmap de presença/ausência\n")
        f.write("- `heatmap_pseudogenization.png` - Heatmap de taxas de pseudogenização\n")
        f.write("- `pca_lysozyme_profiles.png` - PCA dos perfis de lisozimas\n")
        f.write("- `dendrogram_similarity.png` - Dendrograma de similaridade\n")
        f.write("- `boxplots_group_comparison.png` - Comparação entre grupos\n\n")
        
        # SEÇÃO 6: Conclusões
        f.write("## 6. Conclusões\n\n")
        
        # Análise automática de conclusões
        if statistical_tests.get('pseudogene_rate', {}).get('significant', False):
            if path_stats and non_path_stats:
                if path_stats['mean_pseudogene_rate'] > non_path_stats['mean_pseudogene_rate']:
                    f.write("- ✓ **Genomas patogênicos apresentam taxa significativamente MAIOR de pseudogenização de lisozimas.**\n")
                else:
                    f.write("- ✓ **Genomas não-patogênicos apresentam taxa significativamente MAIOR de pseudogenização de lisozimas.**\n")
        else:
            f.write("- Não há diferença estatisticamente significativa na taxa de pseudogenização entre grupos.\n")
        
        if len(protein_categories['core']) > 0:
            f.write(f"- Identificado repertório core de {len(protein_categories['core'])} lisozimas compartilhadas.\n")
        
        if not differential_proteins.empty:
            path_enriched = len(differential_proteins[differential_proteins['enriched_in'] == 'pathogenic'])
            non_path_enriched = len(differential_proteins[differential_proteins['enriched_in'] == 'non-pathogenic'])
            f.write(f"- {path_enriched} lisozimas enriquecidas em patogênicos, {non_path_enriched} em não-patogênicos.\n")
        
        f.write("\n---\n\n")
        f.write("*Relatório gerado automaticamente pelo pipeline de análise comparativa de lisozimas.*\n")
    
    logger.info(f"✓ Relatório executivo salvo: {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Advanced Comparative Lysozyme Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--batch-output',
        type=Path,
        required=True,
        help='Batch pipeline output directory'
    )
    
    parser.add_argument(
        '--metadata',
        type=Path,
        required=True,
        help='TSV file with metadata (genome_id, pathogenicity, ...)'
    )
    
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=None,
        help='Output directory (default: <batch-output>/advanced_comparative)'
    )
    
    parser.add_argument(
        '--alpha',
        type=float,
        default=0.05,
        help='Significance level for statistical tests (default: 0.05)'
    )
    
    args = parser.parse_args()
    
    # Define output dir
    if args.output_dir is None:
        args.output_dir = args.batch_output / "advanced_comparative"
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("=" * 80)
    logger.info("ADVANCED LYSOZYME COMPARATIVE ANALYSIS")
    logger.info("=" * 80)
    
    try:
        # 1. Load data
        logger.info("\n[1/7] Loading data...")
        aggregated_df = load_batch_results(args.batch_output)
        metadata_df = load_metadata(args.metadata)
        
        # 2. Validate metadata
        logger.info("\n[2/7] Validating metadata...")
        is_valid, errors = validate_metadata(metadata_df, aggregated_df)
        
        if not is_valid:
            logger.error("Metadata validation errors:")
            for error in errors:
                logger.error(f"  - {error}")
            sys.exit(1)
        
        logger.info("Metadata validated")
        check_group_balance(metadata_df)
        
        # 3. Protein analysis
        logger.info("\n[3/7] Analyzing protein profiles...")
        protein_profiles = calculate_protein_profiles(aggregated_df)
        protein_categories = classify_protein_categories(protein_profiles, core_threshold=90.0)
        presence_matrix = get_protein_genome_matrix(aggregated_df)
        pseudogenization_matrix = get_pseudogenization_matrix(aggregated_df)
        
        save_protein_profiles(protein_profiles, protein_categories, args.output_dir)
        
        # 4. Pathogenicity analysis
        logger.info("\n[4/7] Comparing pathogenic vs non-pathogenic groups...")
        group_stats = calculate_group_statistics(aggregated_df, metadata_df)
        statistical_tests = perform_statistical_tests(aggregated_df, metadata_df)
        differential_proteins = identify_differential_proteins(aggregated_df, metadata_df)
        
        save_pathogenicity_analysis(
            group_stats,
            statistical_tests,
            differential_proteins,
            args.output_dir
        )
        
        # 5. Visualizations
        logger.info("\n[5/7] Generating visualizations...")
        generate_all_visualizations(
            aggregated_df,
            metadata_df,
            presence_matrix,
            pseudogenization_matrix,
            args.output_dir
        )
        
        # 6. Save matrices
        logger.info("\n[6/7] Saving matrices...")
        presence_matrix.to_csv(args.output_dir / "presence_absence_matrix.tsv", sep='\t')
        pseudogenization_matrix.to_csv(args.output_dir / "pseudogenization_matrix.tsv", sep='\t')
        logger.info("Matrices saved")
        
        # 7. Executive report
        logger.info("\n[7/7] Generating executive report...")
        generate_executive_report(
            aggregated_df,
            metadata_df,
            protein_profiles,
            protein_categories,
            group_stats,
            statistical_tests,
            differential_proteins,
            args.output_dir
        )
        
        logger.info("\n" + "=" * 80)
        logger.info("COMPARATIVE ANALYSIS COMPLETED SUCCESSFULLY")
        logger.info(f"Results saved in: {args.output_dir}")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"\nERROR: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
