"""
Módulo para análise comparativa entre genomas patogênicos e não-patogênicos.
Realiza testes estatísticos e identifica diferenças significativas.
"""

import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple
from scipy import stats


logger = logging.getLogger(__name__)


def validate_metadata(
    metadata_df: pd.DataFrame,
    aggregated_df: pd.DataFrame
) -> Tuple[bool, List[str]]:
    """
    Valida arquivo de metadados contra dados do pipeline.
    
    Args:
        metadata_df: DataFrame com metadados
        aggregated_df: DataFrame com resultados agregados
    
    Returns:
        Tuple (is_valid, list_of_errors)
    """
    errors = []
    
    # Verifica colunas obrigatórias
    required_cols = ['genome_id', 'pathogenicity']
    missing_cols = [col for col in required_cols if col not in metadata_df.columns]
    if missing_cols:
        errors.append(f"Colunas obrigatórias ausentes: {', '.join(missing_cols)}")
    
    # Verifica se genome_ids batem
    metadata_genomes = set(metadata_df['genome_id'].unique())
    data_genomes = set(aggregated_df['genome_id'].unique())
    
    missing_in_metadata = data_genomes - metadata_genomes
    missing_in_data = metadata_genomes - data_genomes
    
    if missing_in_metadata:
        errors.append(f"Genomas sem metadados: {', '.join(missing_in_metadata)}")
    
    if missing_in_data:
        logger.warning(f"Genomas nos metadados mas não nos dados: {', '.join(missing_in_data)}")
    
    # Verifica valores de pathogenicity
    if 'pathogenicity' in metadata_df.columns:
        valid_values = {'pathogenic', 'non-pathogenic', 'unknown'}
        invalid = metadata_df[~metadata_df['pathogenicity'].isin(valid_values)]
        if len(invalid) > 0:
            errors.append(f"Valores inválidos em 'pathogenicity': {invalid['pathogenicity'].unique()}")
    
    is_valid = len(errors) == 0
    
    return is_valid, errors


def check_group_balance(metadata_df: pd.DataFrame, min_group_size: int = 3) -> None:
    """
    Verifica balanceamento de grupos e emite warnings.
    
    Args:
        metadata_df: DataFrame com metadados
        min_group_size: Tamanho mínimo recomendado por grupo
    """
    pathogenic = metadata_df[metadata_df['pathogenicity'] == 'pathogenic']
    non_pathogenic = metadata_df[metadata_df['pathogenicity'] == 'non-pathogenic']
    
    n_pathogenic = len(pathogenic)
    n_non_pathogenic = len(non_pathogenic)
    
    logger.info(f"Grupos: {n_pathogenic} patogênicos, {n_non_pathogenic} não-patogênicos")
    
    if n_pathogenic < min_group_size or n_non_pathogenic < min_group_size:
        logger.warning(
            f"⚠️  ATENÇÃO: Grupos desbalanceados ou pequenos "
            f"(pathogenic={n_pathogenic}, non-pathogenic={n_non_pathogenic}). "
            f"Poder estatístico pode ser insuficiente para detecção de diferenças significativas."
        )
    
    ratio = max(n_pathogenic, n_non_pathogenic) / min(n_pathogenic, n_non_pathogenic) if min(n_pathogenic, n_non_pathogenic) > 0 else float('inf')
    
    if ratio > 3:
        logger.warning(
            f"⚠️  Razão entre grupos muito alta ({ratio:.1f}:1). "
            f"Considere balancear os grupos para melhores resultados estatísticos."
        )


def calculate_group_statistics(
    aggregated_df: pd.DataFrame,
    metadata_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Calcula estatísticas agregadas por grupo (patogênico vs não-patogênico).
    
    Args:
        aggregated_df: DataFrame com resultados
        metadata_df: DataFrame com metadados
    
    Returns:
        DataFrame com estatísticas por grupo
    """
    logger.info("Calculando estatísticas por grupo...")
    
    # Merge com metadados
    merged = aggregated_df.merge(metadata_df[['genome_id', 'pathogenicity']], on='genome_id', how='left')
    
    stats_list = []
    
    for group in ['pathogenic', 'non-pathogenic']:
        group_data = merged[merged['pathogenicity'] == group]
        
        if len(group_data) == 0:
            continue
        
        n_genomes = group_data['genome_id'].nunique()
        n_regions = len(group_data)
        n_pseudogenes = len(group_data[group_data['is_pseudogene'] == True])
        
        stat = {
            'group': group,
            'n_genomes': n_genomes,
            'total_regions': n_regions,
            'functional_genes': n_regions - n_pseudogenes,
            'pseudogenes': n_pseudogenes,
            'mean_pseudogene_rate': (n_pseudogenes / n_regions * 100) if n_regions > 0 else 0,
            'median_score_density': group_data['score_density'].median(),
            'mean_score_density': group_data['score_density'].mean(),
            'std_score_density': group_data['score_density'].std(),
            'total_frameshifts': group_data['frameshifts'].sum(),
            'total_premature_stops': group_data['premature_stop_codons'].sum(),
            'total_missing_start': group_data['missing_start_codon'].sum(),
            'total_missing_stop': group_data['missing_stop_codon'].sum(),
            'mean_frameshifts_per_region': group_data['frameshifts'].mean(),
            'mean_premature_stops_per_region': group_data['premature_stop_codons'].mean()
        }
        
        stats_list.append(stat)
    
    stats_df = pd.DataFrame(stats_list)
    
    logger.info(f"✓ Estatísticas calculadas para {len(stats_df)} grupos")
    
    return stats_df


def perform_statistical_tests(
    aggregated_df: pd.DataFrame,
    metadata_df: pd.DataFrame
) -> Dict[str, Dict]:
    """
    Realiza testes estatísticos comparando grupos.
    
    Args:
        aggregated_df: DataFrame com resultados
        metadata_df: DataFrame com metadados
    
    Returns:
        Dict com resultados dos testes
    """
    logger.info("Realizando testes estatísticos...")
    
    merged = aggregated_df.merge(metadata_df[['genome_id', 'pathogenicity']], on='genome_id', how='left')
    
    pathogenic = merged[merged['pathogenicity'] == 'pathogenic']
    non_pathogenic = merged[merged['pathogenicity'] == 'non-pathogenic']
    
    results = {}
    
    # Teste Mann-Whitney U para taxa de pseudogenes por genoma
    path_pseudogene_rates = pathogenic.groupby('genome_id')['is_pseudogene'].mean() * 100
    non_path_pseudogene_rates = non_pathogenic.groupby('genome_id')['is_pseudogene'].mean() * 100
    
    if len(path_pseudogene_rates) > 0 and len(non_path_pseudogene_rates) > 0:
        stat, pvalue = stats.mannwhitneyu(path_pseudogene_rates, non_path_pseudogene_rates, alternative='two-sided')
        results['pseudogene_rate'] = {
            'test': 'Mann-Whitney U',
            'statistic': stat,
            'p_value': pvalue,
            'significant': pvalue < 0.05,
            'pathogenic_median': path_pseudogene_rates.median(),
            'non_pathogenic_median': non_path_pseudogene_rates.median()
        }
    
    # Teste Mann-Whitney para score density
    if len(pathogenic) > 0 and len(non_pathogenic) > 0:
        stat, pvalue = stats.mannwhitneyu(
            pathogenic['score_density'],
            non_pathogenic['score_density'],
            alternative='two-sided'
        )
        results['score_density'] = {
            'test': 'Mann-Whitney U',
            'statistic': stat,
            'p_value': pvalue,
            'significant': pvalue < 0.05,
            'pathogenic_median': pathogenic['score_density'].median(),
            'non_pathogenic_median': non_pathogenic['score_density'].median()
        }
    
    # Teste para frameshifts
    path_frameshifts = pathogenic.groupby('genome_id')['frameshifts'].sum()
    non_path_frameshifts = non_pathogenic.groupby('genome_id')['frameshifts'].sum()
    
    if len(path_frameshifts) > 0 and len(non_path_frameshifts) > 0:
        stat, pvalue = stats.mannwhitneyu(path_frameshifts, non_path_frameshifts, alternative='two-sided')
        results['frameshifts'] = {
            'test': 'Mann-Whitney U',
            'statistic': stat,
            'p_value': pvalue,
            'significant': pvalue < 0.05,
            'pathogenic_median': path_frameshifts.median(),
            'non_pathogenic_median': non_path_frameshifts.median()
        }
    
    logger.info(f"✓ {len(results)} testes estatísticos completados")
    
    return results


def identify_differential_proteins(
    aggregated_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    min_presence_diff: float = 20.0
) -> pd.DataFrame:
    """
    Identifica proteínas diferencialmente presentes entre grupos.
    
    Args:
        aggregated_df: DataFrame com resultados
        metadata_df: DataFrame com metadados
        min_presence_diff: Diferença mínima de presença para considerar (%)
    
    Returns:
        DataFrame com proteínas diferenciais
    """
    logger.info("Identificando proteínas diferencialmente presentes...")
    
    merged = aggregated_df.merge(metadata_df[['genome_id', 'pathogenicity']], on='genome_id', how='left')
    
    pathogenic_genomes = set(merged[merged['pathogenicity'] == 'pathogenic']['genome_id'].unique())
    non_pathogenic_genomes = set(merged[merged['pathogenicity'] == 'non-pathogenic']['genome_id'].unique())
    
    n_pathogenic = len(pathogenic_genomes)
    n_non_pathogenic = len(non_pathogenic_genomes)
    
    differential_proteins = []
    
    for protein in aggregated_df['best_protein_id'].unique():
        protein_data = merged[merged['best_protein_id'] == protein]
        
        # Presença em cada grupo
        path_with_protein = len(set(protein_data[protein_data['pathogenicity'] == 'pathogenic']['genome_id']))
        non_path_with_protein = len(set(protein_data[protein_data['pathogenicity'] == 'non-pathogenic']['genome_id']))
        
        path_presence = (path_with_protein / n_pathogenic * 100) if n_pathogenic > 0 else 0
        non_path_presence = (non_path_with_protein / n_non_pathogenic * 100) if n_non_pathogenic > 0 else 0
        
        presence_diff = path_presence - non_path_presence
        
        # Pseudogenização em cada grupo
        path_protein_data = protein_data[protein_data['pathogenicity'] == 'pathogenic']
        non_path_protein_data = protein_data[protein_data['pathogenicity'] == 'non-pathogenic']
        
        path_pseudogene_rate = (path_protein_data['is_pseudogene'].sum() / len(path_protein_data) * 100) if len(path_protein_data) > 0 else 0
        non_path_pseudogene_rate = (non_path_protein_data['is_pseudogene'].sum() / len(non_path_protein_data) * 100) if len(non_path_protein_data) > 0 else 0
        
        if abs(presence_diff) >= min_presence_diff or abs(path_pseudogene_rate - non_path_pseudogene_rate) >= 20:
            differential_proteins.append({
                'protein_id': protein,
                'pathogenic_presence': path_presence,
                'non_pathogenic_presence': non_path_presence,
                'presence_difference': presence_diff,
                'pathogenic_pseudogene_rate': path_pseudogene_rate,
                'non_pathogenic_pseudogene_rate': non_path_pseudogene_rate,
                'pseudogene_rate_difference': path_pseudogene_rate - non_path_pseudogene_rate,
                'enriched_in': 'pathogenic' if presence_diff > 0 else 'non-pathogenic'
            })
    
    diff_df = pd.DataFrame(differential_proteins)
    
    if not diff_df.empty:
        diff_df = diff_df.sort_values('presence_difference', key=abs, ascending=False)
    
    logger.info(f"✓ {len(diff_df)} proteínas diferenciais identificadas")
    
    return diff_df


def save_pathogenicity_analysis(
    group_stats: pd.DataFrame,
    statistical_tests: Dict,
    differential_proteins: pd.DataFrame,
    output_dir
) -> None:
    """
    Salva resultados da análise de patogenicidade.
    
    Args:
        group_stats: Estatísticas por grupo
        statistical_tests: Resultados dos testes
        differential_proteins: Proteínas diferenciais
        output_dir: Diretório de saída
    """
    from pathlib import Path
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Salva estatísticas
    stats_path = output_dir / "group_statistics.tsv"
    group_stats.to_csv(stats_path, sep='\t', index=False)
    logger.info(f"✓ Estatísticas de grupos salvas: {stats_path}")
    
    # Salva proteínas diferenciais
    if not differential_proteins.empty:
        diff_path = output_dir / "differential_proteins.tsv"
        differential_proteins.to_csv(diff_path, sep='\t', index=False)
        logger.info(f"✓ Proteínas diferenciais salvas: {diff_path}")
    
    # Relatório textual
    report_path = output_dir / "pathogenicity_differential.txt"
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("ANÁLISE COMPARATIVA: PATOGÊNICOS vs NÃO-PATOGÊNICOS\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("ESTATÍSTICAS POR GRUPO:\n")
        f.write("-" * 80 + "\n")
        for _, row in group_stats.iterrows():
            f.write(f"\n{row['group'].upper()}:\n")
            f.write(f"  Genomas: {row['n_genomes']}\n")
            f.write(f"  Regiões totais: {row['total_regions']}\n")
            f.write(f"  Taxa de pseudogenes: {row['mean_pseudogene_rate']:.2f}%\n")
            f.write(f"  Score density (média±std): {row['mean_score_density']:.2f}±{row['std_score_density']:.2f}\n")
            f.write(f"  Frameshifts totais: {row['total_frameshifts']}\n")
            f.write(f"  Stop codons prematuros: {row['total_premature_stops']}\n")
            f.write(f"  Start codons ausentes: {row['total_missing_start']}\n")
            f.write(f"  Stop codons ausentes: {row['total_missing_stop']}\n")
        
        f.write("\n\nTESTES ESTATÍSTICOS:\n")
        f.write("-" * 80 + "\n")
        for test_name, result in statistical_tests.items():
            f.write(f"\n{test_name.replace('_', ' ').title()}:\n")
            f.write(f"  Teste: {result['test']}\n")
            f.write(f"  p-value: {result['p_value']:.4e}\n")
            f.write(f"  Significativo (α=0.05): {'SIM' if result['significant'] else 'NÃO'}\n")
            if 'pathogenic_median' in result:
                f.write(f"  Mediana patogênicos: {result['pathogenic_median']:.2f}\n")
                f.write(f"  Mediana não-patogênicos: {result['non_pathogenic_median']:.2f}\n")
        
        if not differential_proteins.empty:
            f.write("\n\nPROTEÍNAS DIFERENCIALMENTE PRESENTES:\n")
            f.write("-" * 80 + "\n")
            f.write(f"Total: {len(differential_proteins)} proteínas\n\n")
            
            for _, row in differential_proteins.head(10).iterrows():
                f.write(f"\n{row['protein_id']}:\n")
                f.write(f"  Presença em patogênicos: {row['pathogenic_presence']:.1f}%\n")
                f.write(f"  Presença em não-patogênicos: {row['non_pathogenic_presence']:.1f}%\n")
                f.write(f"  Diferença: {row['presence_difference']:.1f}% (enriquecido em {row['enriched_in']})\n")
                f.write(f"  Pseudogenização em patogênicos: {row['pathogenic_pseudogene_rate']:.1f}%\n")
                f.write(f"  Pseudogenização em não-patogênicos: {row['non_pathogenic_pseudogene_rate']:.1f}%\n")
        
        f.write("\n" + "=" * 80 + "\n")
    
    logger.info(f"✓ Relatório diferencial salvo: {report_path}")
