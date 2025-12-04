"""
Módulo para geração de visualizações comparativas.
Gera heatmaps, PCA, dendrogramas e boxplots para análise visual.
"""

import logging
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Backend não-interativo
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist


logger = logging.getLogger(__name__)

# Configuração de estilo
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


def plot_presence_absence_heatmap(
    matrix: pd.DataFrame,
    metadata_df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Gera heatmap de presença/ausência proteína × genoma.
    
    Args:
        matrix: DataFrame binário (proteínas × genomas)
        metadata_df: DataFrame com metadados (para colorir por patogenicidade)
        output_path: Caminho para salvar figura
    """
    logger.info("Gerando heatmap de presença/ausência...")
    
    # Mapeia pathogenicity
    pathogenicity_map = metadata_df.set_index('genome_id')['pathogenicity'].to_dict()
    
    # Cores para pathogenicity
    genome_colors = []
    for genome in matrix.columns:
        path = pathogenicity_map.get(genome, 'unknown')
        if path == 'pathogenic':
            genome_colors.append('#d62728')  # vermelho
        elif path == 'non-pathogenic':
            genome_colors.append('#2ca02c')  # verde
        else:
            genome_colors.append('#7f7f7f')  # cinza
    
    # Cria figura
    fig, ax = plt.subplots(figsize=(max(12, len(matrix.columns) * 0.4), max(8, len(matrix) * 0.3)))
    
    sns.heatmap(
        matrix,
        cmap='YlOrRd',
        cbar_kws={'label': 'Presença'},
        ax=ax,
        xticklabels=True,
        yticklabels=True,
        linewidths=0.5,
        linecolor='lightgray'
    )
    
    # Adiciona barra de cores para pathogenicity
    for i, color in enumerate(genome_colors):
        ax.add_patch(plt.Rectangle((i, len(matrix)), 1, 0.5, color=color, clip_on=False))
    
    ax.set_xlabel('Genomas', fontsize=12)
    ax.set_ylabel('Proteínas de Referência', fontsize=12)
    ax.set_title('Presença/Ausência de Lisozimas por Genoma', fontsize=14, fontweight='bold')
    
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    logger.info(f"✓ Heatmap salvo: {output_path}")


def plot_pseudogenization_heatmap(
    matrix: pd.DataFrame,
    metadata_df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Gera heatmap de taxa de pseudogenização proteína × genoma.
    
    Args:
        matrix: DataFrame com taxas de pseudogenização (proteínas × genomas)
        metadata_df: DataFrame com metadados
        output_path: Caminho para salvar figura
    """
    logger.info("Gerando heatmap de pseudogenização...")
    
    # Mapeia pathogenicity
    pathogenicity_map = metadata_df.set_index('genome_id')['pathogenicity'].to_dict()
    
    # Cores para pathogenicity
    genome_colors = []
    for genome in matrix.columns:
        path = pathogenicity_map.get(genome, 'unknown')
        if path == 'pathogenic':
            genome_colors.append('#d62728')
        elif path == 'non-pathogenic':
            genome_colors.append('#2ca02c')
        else:
            genome_colors.append('#7f7f7f')
    
    # Cria figura
    fig, ax = plt.subplots(figsize=(max(12, len(matrix.columns) * 0.4), max(8, len(matrix) * 0.3)))
    
    sns.heatmap(
        matrix,
        cmap='RdYlGn_r',
        cbar_kws={'label': 'Taxa de Pseudogenização (%)'},
        ax=ax,
        xticklabels=True,
        yticklabels=True,
        linewidths=0.5,
        linecolor='lightgray',
        vmin=0,
        vmax=100
    )
    
    # Adiciona barra de cores para pathogenicity
    for i, color in enumerate(genome_colors):
        ax.add_patch(plt.Rectangle((i, len(matrix)), 1, 0.5, color=color, clip_on=False))
    
    ax.set_xlabel('Genomas', fontsize=12)
    ax.set_ylabel('Proteínas de Referência', fontsize=12)
    ax.set_title('Taxa de Pseudogenização de Lisozimas por Genoma', fontsize=14, fontweight='bold')
    
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    logger.info(f"✓ Heatmap de pseudogenização salvo: {output_path}")


def plot_pca(
    aggregated_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Gera PCA dos perfis de lisozimas colorido por patogenicidade.
    
    Args:
        aggregated_df: DataFrame agregado
        metadata_df: DataFrame com metadados
        output_path: Caminho para salvar figura
    """
    logger.info("Gerando PCA...")
    
    # Cria matriz de features por genoma
    features = []
    genome_ids = []
    
    for genome_id in aggregated_df['genome_id'].unique():
        genome_data = aggregated_df[aggregated_df['genome_id'] == genome_id]
        
        feature_vector = [
            len(genome_data),  # número de regiões
            genome_data['is_pseudogene'].mean() * 100,  # taxa de pseudogenes
            genome_data['score_density'].mean(),  # score density médio
            genome_data['frameshifts'].sum(),  # frameshifts totais
            genome_data['premature_stop_codons'].sum(),  # stops prematuros
            genome_data['missing_start_codon'].sum(),  # starts ausentes
            genome_data['missing_stop_codon'].sum(),  # stops ausentes
            genome_data['best_protein_id'].nunique()  # diversidade de proteínas
        ]
        
        features.append(feature_vector)
        genome_ids.append(genome_id)
    
    # Normaliza features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)
    
    # PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(features_scaled)
    
    # Mapeia pathogenicity
    pathogenicity_map = metadata_df.set_index('genome_id')['pathogenicity'].to_dict()
    colors = [pathogenicity_map.get(g, 'unknown') for g in genome_ids]
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    color_map = {'pathogenic': '#d62728', 'non-pathogenic': '#2ca02c', 'unknown': '#7f7f7f'}
    
    for pathogenicity, color in color_map.items():
        mask = [c == pathogenicity for c in colors]
        if any(mask):
            ax.scatter(
                pca_result[mask, 0],
                pca_result[mask, 1],
                c=color,
                label=pathogenicity.replace('_', ' ').title(),
                s=100,
                alpha=0.7,
                edgecolors='black'
            )
    
    # Anota genomas
    for i, genome_id in enumerate(genome_ids):
        ax.annotate(
            genome_id,
            (pca_result[i, 0], pca_result[i, 1]),
            fontsize=8,
            alpha=0.7,
            xytext=(5, 5),
            textcoords='offset points'
        )
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
    ax.set_title('PCA - Perfis de Lisozimas', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    logger.info(f"✓ PCA salvo: {output_path}")


def plot_dendrogram(
    presence_matrix: pd.DataFrame,
    metadata_df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Gera dendrograma hierárquico de similaridade entre genomas.
    
    Args:
        presence_matrix: Matriz de presença/ausência
        metadata_df: DataFrame com metadados
        output_path: Caminho para salvar figura
    """
    logger.info("Gerando dendrograma...")
    
    # Transpõe para ter genomas como linhas
    genome_matrix = presence_matrix.T
    
    # Calcula linkage
    linkage_matrix = linkage(genome_matrix, method='ward', metric='euclidean')
    
    # Mapeia pathogenicity para cores
    pathogenicity_map = metadata_df.set_index('genome_id')['pathogenicity'].to_dict()
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    dendro = dendrogram(
        linkage_matrix,
        labels=genome_matrix.index.tolist(),
        ax=ax,
        orientation='right'
    )
    
    # Colore labels por pathogenicity
    color_map = {'pathogenic': '#d62728', 'non-pathogenic': '#2ca02c', 'unknown': '#7f7f7f'}
    
    ax.set_xlabel('Distância', fontsize=12)
    ax.set_title('Dendrograma - Similaridade de Perfis de Lisozimas', fontsize=14, fontweight='bold')
    
    # Colore os labels
    ylbls = ax.get_ymajorticklabels()
    for lbl in ylbls:
        genome_id = lbl.get_text()
        path = pathogenicity_map.get(genome_id, 'unknown')
        lbl.set_color(color_map.get(path, '#7f7f7f'))
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    logger.info(f"✓ Dendrograma salvo: {output_path}")


def plot_group_comparisons(
    aggregated_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Gera boxplots comparando grupos patogênicos vs não-patogênicos.
    
    Args:
        aggregated_df: DataFrame agregado
        metadata_df: DataFrame com metadados
        output_path: Caminho para salvar figura
    """
    logger.info("Gerando boxplots comparativos...")
    
    # Merge com metadados
    merged = aggregated_df.merge(metadata_df[['genome_id', 'pathogenicity']], on='genome_id', how='left')
    merged = merged[merged['pathogenicity'].isin(['pathogenic', 'non-pathogenic'])]
    
    # Calcula métricas por genoma
    genome_metrics = merged.groupby(['genome_id', 'pathogenicity']).agg({
        'is_pseudogene': lambda x: (x.sum() / len(x)) * 100,
        'score_density': 'mean',
        'frameshifts': 'sum',
        'premature_stop_codons': 'sum'
    }).reset_index()
    
    genome_metrics.columns = ['genome_id', 'pathogenicity', 'pseudogene_rate', 'mean_score_density', 'total_frameshifts', 'total_premature_stops']
    
    # Cria subplot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Taxa de pseudogenes
    sns.boxplot(
        data=genome_metrics,
        x='pathogenicity',
        y='pseudogene_rate',
        palette={'pathogenic': '#d62728', 'non-pathogenic': '#2ca02c'},
        ax=axes[0, 0]
    )
    axes[0, 0].set_title('Taxa de Pseudogenes', fontweight='bold')
    axes[0, 0].set_ylabel('Pseudogenes (%)')
    axes[0, 0].set_xlabel('')
    
    # 2. Score density
    sns.boxplot(
        data=genome_metrics,
        x='pathogenicity',
        y='mean_score_density',
        palette={'pathogenic': '#d62728', 'non-pathogenic': '#2ca02c'},
        ax=axes[0, 1]
    )
    axes[0, 1].set_title('Score Density Médio', fontweight='bold')
    axes[0, 1].set_ylabel('Score Density')
    axes[0, 1].set_xlabel('')
    
    # 3. Frameshifts
    sns.boxplot(
        data=genome_metrics,
        x='pathogenicity',
        y='total_frameshifts',
        palette={'pathogenic': '#d62728', 'non-pathogenic': '#2ca02c'},
        ax=axes[1, 0]
    )
    axes[1, 0].set_title('Frameshifts Totais', fontweight='bold')
    axes[1, 0].set_ylabel('Número de Frameshifts')
    axes[1, 0].set_xlabel('Grupo')
    
    # 4. Stop codons prematuros
    sns.boxplot(
        data=genome_metrics,
        x='pathogenicity',
        y='total_premature_stops',
        palette={'pathogenic': '#d62728', 'non-pathogenic': '#2ca02c'},
        ax=axes[1, 1]
    )
    axes[1, 1].set_title('Stop Codons Prematuros', fontweight='bold')
    axes[1, 1].set_ylabel('Número de Stops Prematuros')
    axes[1, 1].set_xlabel('Grupo')
    
    plt.suptitle('Comparação: Patogênicos vs Não-Patogênicos', fontsize=16, fontweight='bold', y=1.00)
    plt.tight_layout()
    
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    logger.info(f"✓ Boxplots salvos: {output_path}")


def generate_all_visualizations(
    aggregated_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    presence_matrix: pd.DataFrame,
    pseudogenization_matrix: pd.DataFrame,
    output_dir: Path
) -> None:
    """
    Gera todas as visualizações.
    
    Args:
        aggregated_df: DataFrame agregado
        metadata_df: DataFrame com metadados
        presence_matrix: Matriz de presença/ausência
        pseudogenization_matrix: Matriz de pseudogenização
        output_dir: Diretório para salvar figuras
    """
    logger.info("Gerando todas as visualizações...")
    
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Heatmap presença/ausência
    plot_presence_absence_heatmap(
        presence_matrix,
        metadata_df,
        figures_dir / "heatmap_presence_absence.png"
    )
    
    # 2. Heatmap pseudogenização
    plot_pseudogenization_heatmap(
        pseudogenization_matrix,
        metadata_df,
        figures_dir / "heatmap_pseudogenization.png"
    )
    
    # 3. PCA
    plot_pca(
        aggregated_df,
        metadata_df,
        figures_dir / "pca_lysozyme_profiles.png"
    )
    
    # 4. Dendrograma
    plot_dendrogram(
        presence_matrix,
        metadata_df,
        figures_dir / "dendrogram_similarity.png"
    )
    
    # 5. Boxplots comparativos
    plot_group_comparisons(
        aggregated_df,
        metadata_df,
        figures_dir / "boxplots_group_comparison.png"
    )
    
    logger.info(f"✓ Todas as visualizações salvas em: {figures_dir}")
