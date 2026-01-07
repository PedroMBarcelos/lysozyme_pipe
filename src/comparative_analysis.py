"""
Module for comparative analysis of multiple genomes.
Aggregates results from different genomes and generates consolidated statistics.
"""

import logging
import pandas as pd
from pathlib import Path
from typing import List, Dict
from collections import defaultdict

from src.pseudogene_detection import PseudogeneAnnotation


logger = logging.getLogger(__name__)


def aggregate_all_genomes(results_list: List[pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenate results from all genomes into a single DataFrame.
    
    Args:
        results_list: List of DataFrames, each representing a genome
    
    Returns:
        Aggregated DataFrame with all genomes
    """
    if not results_list:
        logger.warning("No results to aggregate")
        return pd.DataFrame()
    
    logger.debug(f"Aggregating results from {len(results_list)} genomes...")
    
    aggregated_df = pd.concat(results_list, ignore_index=True)
    
    logger.info(f"Aggregation complete. Total regions: {len(aggregated_df)}")
    
    return aggregated_df


def generate_summary_stats(aggregated_df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate summary statistics per genome.
    
    Args:
        aggregated_df: DataFrame with aggregated results
    
    Returns:
        DataFrame with statistics per genome
    """
    logger.debug("Generating summary statistics per genome...")
    
    if aggregated_df.empty:
        return pd.DataFrame()
    
    stats = []
    
    for genome_id in aggregated_df['genome_id'].unique():
        genome_data = aggregated_df[aggregated_df['genome_id'] == genome_id]
        
        stat = {
            'genome_id': genome_id,
            'total_regions': len(genome_data),
            'functional_genes': len(genome_data[genome_data['is_pseudogene'] == False]),
            'pseudogenes': len(genome_data[genome_data['is_pseudogene'] == True]),
            'functional_possible': len(genome_data[(genome_data['is_pseudogene'] == False) & (genome_data['is_small_orf'] == False)]),
            'functional_small': len(genome_data[(genome_data['is_pseudogene'] == False) & (genome_data['is_small_orf'] == True)]),
            'pseudogene_detected': len(genome_data[(genome_data['is_pseudogene'] == True) & (genome_data['is_small_orf'] == False)]),
            'pseudogene_small': len(genome_data[(genome_data['is_pseudogene'] == True) & (genome_data['is_small_orf'] == True)]),
            'pseudogene_rate': len(genome_data[genome_data['is_pseudogene'] == True]) / len(genome_data) * 100 if len(genome_data) > 0 else 0,
            'total_stop_codons': genome_data['premature_stop_codons'].sum(),
            'total_frameshifts': genome_data['frameshifts'].sum(),
            'total_missing_start': genome_data['missing_start_codon'].sum(),
            'total_missing_stop': genome_data['missing_stop_codon'].sum(),
            'total_substitutions': genome_data['non_synonymous_substitutions'].sum(),
            'total_indels': genome_data['in_frame_indels'].sum(),
            'mean_score_density': genome_data['score_density'].mean(),
            'unique_proteins': genome_data['best_protein_id'].nunique()
        }
        
        stats.append(stat)
    
    stats_df = pd.DataFrame(stats)
    
    logger.debug(f"Statistics generated for {len(stats)} genomes")
    
    return stats_df


def create_presence_absence_matrix(aggregated_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create presence/absence matrix of lysozymes per genome.
    
    Rows: Reference proteins (lysozymes)
    Columns: Genomes
    Values: 1 (present), 0 (absent)
    
    Args:
        aggregated_df: DataFrame with aggregated results
    
    Returns:
        DataFrame with presence/absence matrix
    """
    logger.debug("Creating presence/absence matrix...")
    
    if aggregated_df.empty:
        return pd.DataFrame()
    
    # Create empty matrix
    genomes = sorted(aggregated_df['genome_id'].unique())
    proteins = sorted(aggregated_df['best_protein_id'].unique())
    
    matrix = pd.DataFrame(0, index=proteins, columns=genomes)
    
    # Fill matrix
    for _, row in aggregated_df.iterrows():
        genome = row['genome_id']
        protein = row['best_protein_id']
        matrix.at[protein, genome] = 1
    
    logger.debug(f"Matrix created: {len(proteins)} proteins × {len(genomes)} genomes")
    
    return matrix


def generate_comparative_report(
    aggregated_df: pd.DataFrame,
    output_dir: Path
) -> None:
    """
    Generate complete comparative report with tables and statistics.
    
    Args:
        aggregated_df: DataFrame with aggregated results
        output_dir: Directory to save files
    """
    logger.debug("Generating comparative report...")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Save complete aggregated table
    aggregated_path = output_dir / "all_genomes_aggregated.tsv"
    aggregated_df.to_csv(aggregated_path, sep='\t', index=False)
    logger.debug(f"Aggregated table saved: {aggregated_path}")
    
    # 2. Generate and save summary statistics
    summary_stats = generate_summary_stats(aggregated_df)
    summary_path = output_dir / "summary_stats.csv"
    summary_stats.to_csv(summary_path, index=False)
    logger.debug(f"Summary statistics saved: {summary_path}")
    
    # 3. Generate and save presence/absence matrix
    presence_absence = create_presence_absence_matrix(aggregated_df)
    matrix_path = output_dir / "presence_absence_matrix.tsv"
    presence_absence.to_csv(matrix_path, sep='\t')
    logger.debug(f"Presence/absence matrix saved: {matrix_path}")
    
    # 4. Generate text summary report
    report_path = output_dir / "comparative_summary.txt"
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("MULTI-GENOME COMPARATIVE REPORT\n")
        f.write("Lysozyme Annotation Pipeline\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total genomes analyzed: {len(aggregated_df['genome_id'].unique())}\n")
        f.write(f"Total regions identified: {len(aggregated_df)}\n")
        f.write(f"Total reference proteins: {aggregated_df['best_protein_id'].nunique()}\n\n")
        
        f.write("STATISTICS PER GENOME:\n")
        f.write("-" * 80 + "\n")
        
        for _, row in summary_stats.iterrows():
            f.write(f"\n{row['genome_id']}:\n")
            f.write(f"  Total regions: {row['total_regions']}\n")
            f.write(f"  Functional genes: {row['functional_genes']} ({100-row['pseudogene_rate']:.1f}%)\n")
            f.write(f"    - Possible Genes: {row['functional_possible']}\n")
            f.write(f"    - Small ORFs: {row['functional_small']}\n")
            f.write(f"  Pseudogenes: {row['pseudogenes']} ({row['pseudogene_rate']:.1f}%)\n")
            f.write(f"    - Detected Pseudogenes: {row['pseudogene_detected']}\n")
            f.write(f"    - Small ORFs: {row['pseudogene_small']}\n")
            f.write(f"  Premature stop codons: {row['total_stop_codons']}\n")
            f.write(f"  Frameshifts: {row['total_frameshifts']}\n")
            f.write(f"  Missing start codons: {row['total_missing_start']}\n")
            f.write(f"  Missing stop codons: {row['total_missing_stop']}\n")
            f.write(f"  Mean score density: {row['mean_score_density']:.2f}\n")
            f.write(f"  Unique proteins detected: {row['unique_proteins']}\n")
        
        f.write("\n" + "=" * 80 + "\n")
    
    logger.debug(f"Text report saved: {report_path}")
    
    # 5. Generate visualizations
    try:
        from src.visualization import generate_all_visualizations
        
        # Check if metadata file exists (look in multiple locations)
        metadata_paths = [
            output_dir.parent / "metadata.tsv",  # Same parent as comparative_analysis
            Path("downloaded_genomes") / "metadata.tsv",  # Default download location
            Path("genomes") / "metadata.tsv",  # Alternative location
        ]
        
        metadata_df = None
        for metadata_path in metadata_paths:
            if metadata_path.exists():
                logger.info(f"Found metadata file: {metadata_path}")
                metadata_df = pd.read_csv(metadata_path, sep='\t', dtype={'genome_id': str})
                break
        
        if metadata_df is not None:
            
            # Create presence/absence matrix for visualizations
            presence_matrix = create_presence_absence_matrix(aggregated_df)
            
            # Create pseudogenization matrix
            pseudogenization_matrix = presence_matrix.copy()
            for genome_id in aggregated_df['genome_id'].unique():
                genome_data = aggregated_df[aggregated_df['genome_id'] == genome_id]
                for protein_id in presence_matrix.index:
                    protein_data = genome_data[genome_data['best_protein_id'] == protein_id]
                    if not protein_data.empty:
                        # Mark as 2 if pseudogene, 1 if functional
                        is_pseudo = protein_data.iloc[0]['is_pseudogene']
                        pseudogenization_matrix.loc[protein_id, genome_id] = 2 if is_pseudo else 1
            
            generate_all_visualizations(
                aggregated_df,
                metadata_df,
                presence_matrix,
                pseudogenization_matrix,
                output_dir
            )
        else:
            logger.warning(f"Metadata file not found at {metadata_path}, skipping visualizations")
    except ImportError as e:
        logger.warning(f"Visualization dependencies not available: {e}")
    except Exception as e:
        logger.error(f"Error generating visualizations: {e}", exc_info=True)
    
    logger.info("Complete comparative report generated")


def compare_pseudogene_rates(summary_stats: pd.DataFrame) -> Dict[str, float]:
    """
    Compare pseudogenization rates across genomes.
    
    Args:
        summary_stats: DataFrame with statistics per genome
    
    Returns:
        Dictionary with comparative statistics
    """
    if summary_stats.empty:
        return {}
    
    return {
        'mean_pseudogene_rate': summary_stats['pseudogene_rate'].mean(),
        'min_pseudogene_rate': summary_stats['pseudogene_rate'].min(),
        'max_pseudogene_rate': summary_stats['pseudogene_rate'].max(),
        'std_pseudogene_rate': summary_stats['pseudogene_rate'].std(),
        'genome_with_most_pseudogenes': summary_stats.loc[summary_stats['pseudogenes'].idxmax(), 'genome_id'],
        'genome_with_least_pseudogenes': summary_stats.loc[summary_stats['pseudogenes'].idxmin(), 'genome_id']
    }
