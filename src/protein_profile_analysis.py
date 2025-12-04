"""
Módulo para análise de perfis de proteínas de referência.
Analisa a distribuição e características de cada lisozima de referência nos genomas.
"""

import logging
import pandas as pd
from typing import Dict, List, Tuple
from collections import defaultdict


logger = logging.getLogger(__name__)


def calculate_protein_profiles(aggregated_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calcula perfis de distribuição para cada proteína de referência.
    
    Args:
        aggregated_df: DataFrame com resultados agregados de todos os genomas
    
    Returns:
        DataFrame com estatísticas por proteína
    """
    logger.info("Calculando perfis de proteínas de referência...")
    
    if aggregated_df.empty:
        logger.warning("DataFrame vazio, sem dados para análise")
        return pd.DataFrame()
    
    total_genomes = aggregated_df['genome_id'].nunique()
    profiles = []
    
    for protein_id in aggregated_df['best_protein_id'].unique():
        protein_data = aggregated_df[aggregated_df['best_protein_id'] == protein_id]
        
        # Genomas com esta proteína
        genomes_with_protein = protein_data['genome_id'].nunique()
        presence_rate = (genomes_with_protein / total_genomes) * 100
        
        # Taxa de pseudogenização
        total_regions = len(protein_data)
        pseudogenes = len(protein_data[protein_data['is_pseudogene'] == True])
        pseudogene_rate = (pseudogenes / total_regions) * 100 if total_regions > 0 else 0
        
        # Estatísticas de score
        mean_score_density = protein_data['score_density'].mean()
        median_score_density = protein_data['score_density'].median()
        std_score_density = protein_data['score_density'].std()
        
        # Tipos de mutações
        total_frameshifts = protein_data['frameshifts'].sum()
        total_premature_stops = protein_data['premature_stop_codons'].sum()
        total_missing_start = protein_data['missing_start_codon'].sum()
        total_missing_stop = protein_data['missing_stop_codon'].sum()
        
        profile = {
            'protein_id': protein_id,
            'genomes_with_protein': genomes_with_protein,
            'total_genomes': total_genomes,
            'presence_rate': presence_rate,
            'total_regions': total_regions,
            'functional_regions': total_regions - pseudogenes,
            'pseudogene_regions': pseudogenes,
            'pseudogene_rate': pseudogene_rate,
            'mean_score_density': mean_score_density,
            'median_score_density': median_score_density,
            'std_score_density': std_score_density,
            'total_frameshifts': total_frameshifts,
            'total_premature_stops': total_premature_stops,
            'total_missing_start': total_missing_start,
            'total_missing_stop': total_missing_stop
        }
        
        profiles.append(profile)
    
    profiles_df = pd.DataFrame(profiles)
    profiles_df = profiles_df.sort_values('presence_rate', ascending=False)
    
    logger.info(f"✓ Perfis calculados para {len(profiles_df)} proteínas")
    
    return profiles_df


def classify_protein_categories(
    profiles_df: pd.DataFrame,
    core_threshold: float = 90.0
) -> Dict[str, List[str]]:
    """
    Classifica proteínas em categorias baseado em taxa de presença.
    
    Args:
        profiles_df: DataFrame com perfis de proteínas
        core_threshold: Threshold para considerar proteína "core" (padrão: 90%)
    
    Returns:
        Dict com listas de protein_ids por categoria
    """
    logger.info(f"Classificando proteínas (core threshold={core_threshold}%)...")
    
    categories = {
        'core': [],        # ≥90% dos genomas
        'accessory': [],   # <90% dos genomas, >1 genoma
        'unique': []       # apenas 1 genoma
    }
    
    for _, row in profiles_df.iterrows():
        protein_id = row['protein_id']
        presence_rate = row['presence_rate']
        genomes_with_protein = row['genomes_with_protein']
        
        if presence_rate >= core_threshold:
            categories['core'].append(protein_id)
        elif genomes_with_protein > 1:
            categories['accessory'].append(protein_id)
        else:
            categories['unique'].append(protein_id)
    
    logger.info(f"  Core: {len(categories['core'])} proteínas")
    logger.info(f"  Accessory: {len(categories['accessory'])} proteínas")
    logger.info(f"  Unique: {len(categories['unique'])} proteínas")
    
    return categories


def get_protein_genome_matrix(aggregated_df: pd.DataFrame) -> pd.DataFrame:
    """
    Cria matriz de presença/ausência proteína × genoma.
    
    Args:
        aggregated_df: DataFrame agregado
    
    Returns:
        DataFrame binário (1=presente, 0=ausente)
    """
    logger.info("Criando matriz proteína × genoma...")
    
    genomes = sorted(aggregated_df['genome_id'].unique())
    proteins = sorted(aggregated_df['best_protein_id'].unique())
    
    matrix = pd.DataFrame(0, index=proteins, columns=genomes)
    
    for _, row in aggregated_df.iterrows():
        protein = row['best_protein_id']
        genome = row['genome_id']
        matrix.at[protein, genome] = 1
    
    logger.info(f"✓ Matriz criada: {len(proteins)} proteínas × {len(genomes)} genomas")
    
    return matrix


def get_pseudogenization_matrix(aggregated_df: pd.DataFrame) -> pd.DataFrame:
    """
    Cria matriz de pseudogenização proteína × genoma.
    
    Args:
        aggregated_df: DataFrame agregado
    
    Returns:
        DataFrame com taxa de pseudogenização (0-100%) ou NaN se ausente
    """
    logger.info("Criando matriz de pseudogenização...")
    
    genomes = sorted(aggregated_df['genome_id'].unique())
    proteins = sorted(aggregated_df['best_protein_id'].unique())
    
    matrix = pd.DataFrame(float('nan'), index=proteins, columns=genomes)
    
    for protein in proteins:
        for genome in genomes:
            data = aggregated_df[
                (aggregated_df['best_protein_id'] == protein) & 
                (aggregated_df['genome_id'] == genome)
            ]
            
            if len(data) > 0:
                pseudogene_rate = (data['is_pseudogene'].sum() / len(data)) * 100
                matrix.at[protein, genome] = pseudogene_rate
    
    logger.info(f"✓ Matriz de pseudogenização criada")
    
    return matrix


def save_protein_profiles(
    profiles_df: pd.DataFrame,
    categories: Dict[str, List[str]],
    output_dir
) -> None:
    """
    Salva perfis de proteínas e categorias.
    
    Args:
        profiles_df: DataFrame com perfis
        categories: Dict com categorias de proteínas
        output_dir: Diretório de saída
    """
    from pathlib import Path
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Salva perfis completos
    profiles_path = output_dir / "protein_profiles.tsv"
    profiles_df.to_csv(profiles_path, sep='\t', index=False)
    logger.info(f"✓ Perfis salvos: {profiles_path}")
    
    # Salva categorias
    categories_path = output_dir / "protein_categories.txt"
    with open(categories_path, 'w') as f:
        f.write("CATEGORIZAÇÃO DE PROTEÍNAS DE REFERÊNCIA\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"CORE LYSOZYMES (presentes em ≥90% dos genomas):\n")
        f.write(f"  Total: {len(categories['core'])} proteínas\n")
        for protein in categories['core']:
            f.write(f"  - {protein}\n")
        
        f.write(f"\nACCESSORY LYSOZYMES (presentes em múltiplos genomas, <90%):\n")
        f.write(f"  Total: {len(categories['accessory'])} proteínas\n")
        for protein in categories['accessory']:
            f.write(f"  - {protein}\n")
        
        f.write(f"\nUNIQUE LYSOZYMES (presentes em apenas 1 genoma):\n")
        f.write(f"  Total: {len(categories['unique'])} proteínas\n")
        for protein in categories['unique']:
            f.write(f"  - {protein}\n")
    
    logger.info(f"✓ Categorias salvas: {categories_path}")
