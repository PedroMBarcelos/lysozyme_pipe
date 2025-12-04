#!/usr/bin/env python3
"""
Main lysozyme annotation and pseudogene detection pipeline.
Orchestrates all pipeline steps in an integrated manner.
"""

import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime

# Importa todos os módulos do pipeline
from src.config import (
    DEFAULT_BLAST_PARAMS,
    DEFAULT_SSEARCH_PARAMS,
    DEFAULT_BEDTOOLS_PARAMS,
    MIN_IDENTITY_THRESHOLD,
    MIN_BLOSUM62_SCORE
)
from src.dependencies import verify_and_install_dependencies
from src.blast_search import run_blast_pipeline_step
from src.blast_filter import run_filtering_step
from src.ssearch_realign import realign_filtered_hits
from src.bedtools_merge import merge_blast_hits, save_merged_regions
from src.score_density import annotate_regions_with_best_proteins, save_region_annotations
from src.pseudogene_detection import (
    annotate_pseudogenes,
    save_pseudogene_annotations,
    generate_summary_report
)


def setup_logging(log_file: Path = None, verbose: bool = False) -> None:
    """
    Configure logging system.
    
    Args:
        log_file: File to save logs (optional)
        verbose: If True, displays DEBUG messages
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    
    # Message format
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    
    # Basic configuration
    handlers = [logging.StreamHandler(sys.stdout)]
    
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=log_level,
        format=log_format,
        datefmt=date_format,
        handlers=handlers
    )


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(
        description='Pipeline de anotação de lisozimas e detecção de pseudogenes em genomas de E. coli',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:
  
  # Execução básica
  python pipeline.py -g genome.fasta -l lysozymes.fasta -o results/
  
  # Com parâmetros customizados
  python pipeline.py -g genome.fasta -l lysozymes.fasta -o results/ \\
                     --min-identity 25 --min-score 150 -v
  
  # Apenas etapas específicas
  python pipeline.py -g genome.fasta -l lysozymes.fasta -o results/ \\
                     --skip-ssearch
        """
    )
    
    # Argumentos obrigatórios
    required = parser.add_argument_group('argumentos obrigatórios')
    
    # Modo single genome OU modo lote
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-g', '--genome',
        type=Path,
        metavar='FASTA',
        help='FASTA file of a single E. coli genome'
    )
    input_group.add_argument(
        '--input-dir',
        type=Path,
        metavar='DIR',
        help='Directory with multiple FASTA files (batch mode)'
    )
    
    required.add_argument(
        '-l', '--lysozymes',
        type=Path,
        required=True,
        metavar='FASTA',
        help='FASTA file with reference lysozymes (UniProtKB/Swiss-Prot)'
    )
    required.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        metavar='DIR',
        help='Output directory for results'
    )
    
    # Filtering parameters
    filtering = parser.add_argument_group('filtering parameters')
    filtering.add_argument(
        '--min-identity',
        type=float,
        default=MIN_IDENTITY_THRESHOLD,
        metavar='FLOAT',
        help=f'Minimum percent identity (default: {MIN_IDENTITY_THRESHOLD})'
    )
    filtering.add_argument(
        '--min-score',
        type=int,
        default=MIN_BLOSUM62_SCORE,
        metavar='INT',
        help=f'Minimum BLOSUM62 score (default: {MIN_BLOSUM62_SCORE})'
    )
    filtering.add_argument(
        '--min-disablements',
        type=int,
        default=1,
        metavar='INT',
        help='Minimum number of mutations to classify as pseudogene (default: 1)'
    )
    
    # Execution options
    execution = parser.add_argument_group('execution options')
    execution.add_argument(
        '--num-threads',
        type=int,
        default=None,
        metavar='INT',
        help='Número de threads para SSEARCH paralelo (padrão: CPU_COUNT-1)'
    )
    execution.add_argument(
        '--skip-ssearch',
        action='store_true',
        help='Pula a etapa de realinhamento com SSEARCH (mais rápido, menos preciso)'
    )
    execution.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Exibe mensagens de debug detalhadas'
    )
    execution.add_argument(
        '--log-file',
        type=Path,
        metavar='FILE',
        help='Arquivo para salvar logs (opcional)'
    )
    
    return parser.parse_args()


def validate_inputs(genome_path: Path, lysozymes_path: Path) -> None:
    """
    Valida os arquivos de entrada.
    
    Args:
        genome_path: Caminho para o arquivo do genoma
        lysozymes_path: Caminho para o arquivo de lisozimas
    
    Raises:
        FileNotFoundError: Se algum arquivo não existir
    """
    if not genome_path.exists():
        raise FileNotFoundError(f"Arquivo de genoma não encontrado: {genome_path}")
    
    if not lysozymes_path.exists():
        raise FileNotFoundError(f"Arquivo de lisozimas não encontrado: {lysozymes_path}")


def run_pipeline(
    genome_fasta: Path,
    lysozyme_fasta: Path,
    output_dir: Path,
    genome_id: str = "unknown",
    min_identity: float = MIN_IDENTITY_THRESHOLD,
    min_score: int = MIN_BLOSUM62_SCORE,
    min_disablements: int = 1,
    skip_ssearch: bool = False,
    num_threads: int = None,
    deps: dict = None
) -> None:
    """
    Executa o pipeline completo de anotação de lisozimas.
    
    Args:
        genome_fasta: Arquivo FASTA do genoma
        lysozyme_fasta: Arquivo FASTA com lisozimas de referência
        output_dir: Output directory
        genome_id: Genome identifier
        min_identity: Minimum identity for filtering
        min_score: Minimum score for filtering
        min_disablements: Minimum mutations to classify as pseudogene
        skip_ssearch: If True, skip SSEARCH realignment
        num_threads: Number of threads for parallelization
        deps: Dependencies dictionary from verify_and_install_dependencies()
    """
    logger = logging.getLogger(__name__)
    
    # If deps not provided, verify dependencies (for standalone usage)
    if deps is None:
        deps = verify_and_install_dependencies()
    
    # Configure BLAST threads
    from src.config import DEFAULT_BLAST_PARAMS
    blast_threads = num_threads if num_threads is not None else deps['num_threads']
    DEFAULT_BLAST_PARAMS.num_threads = blast_threads
    # Create output directories
    output_dir.mkdir(parents=True, exist_ok=True)
    blast_dir = output_dir / "blast"
    merge_dir = output_dir / "merge"
    final_dir = output_dir / "final"
    
    # STEP 1: BLAST search
    logger.info(f"[1/3] BLAST search ({blast_threads} threads)...")
    
    blast_output = run_blast_pipeline_step(
        genome_fasta,
        lysozyme_fasta,
        blast_dir,
        DEFAULT_BLAST_PARAMS,
        deps['makeblastdb'],
        deps['tblastn']
    )
    
    # STEP 2: Quality filtering (silent)
    filtered_output = blast_dir / "filtered_hits.tsv"
    filtered_hits = run_filtering_step(
        blast_output,
        filtered_output,
        min_identity,
        min_score,
        genome_id=genome_id
    )
    
    if not filtered_hits:
        logger.error("No hits passed filtering")
        return
    
    # STEP 3: SSEARCH realignment (silent)
    ssearch_dir = output_dir / "ssearch"
    ssearch_dir.mkdir(parents=True, exist_ok=True)
    
    from src.ssearch_realign import realign_filtered_hits_parallel
    
    ssearch_alignments = realign_filtered_hits_parallel(
        filtered_hits=filtered_hits,
        query_fasta_path=lysozyme_fasta,
        genome_fasta_path=genome_fasta,
        output_dir=ssearch_dir,
        ssearch_path=deps['ssearch36'],
        num_threads=num_threads
    )
    
    logger.debug(f"SSEARCH realignments complete: {len(ssearch_alignments)}")
    
    # Filter realignments by E-value
    evalue_threshold = 1e-7
    filtered_ssearch = {key: aln for key, aln in ssearch_alignments.items() 
                       if aln.evalue <= evalue_threshold}
    
    logger.debug(f"Filtering SSEARCH by E-value <= {evalue_threshold:.0e}")
    logger.debug(f"  - Before filtering: {len(ssearch_alignments)} realignments")
    logger.debug(f"  - After filtering: {len(filtered_ssearch)} realignments")
    
    # Save SSEARCH results
    ssearch_output = ssearch_dir / "ssearch_realignments.tsv"
    with open(ssearch_output, 'w') as f:
        f.write("hit_key\tquery_id\tsubject_id\tidentity\tevalue\tbit_score\n")
        for key, aln in filtered_ssearch.items():
            f.write(f"{key}\t{aln.query_id}\t{aln.subject_id}\t"
                   f"{aln.identity:.2f}\t{aln.evalue:.2e}\t{aln.bit_score:.2f}\n")
    
    # Update filtered hits to use only those that passed SSEARCH
    ssearch_hit_keys = set(filtered_ssearch.keys())
    filtered_hits_after_ssearch = []
    for hit in filtered_hits:
        hit_key = f"{hit.qseqid}_{hit.sseqid}_{hit.sstart}_{hit.send}"
        if hit_key in ssearch_hit_keys:
            filtered_hits_after_ssearch.append(hit)
    
    if not filtered_hits_after_ssearch:
        logger.error("No hits passed SSEARCH filtering")
        return
    
    filtered_hits = filtered_hits_after_ssearch
    
    # STEP 2: Region merging
    logger.info(f"[2/3] Merging {len(filtered_hits)} regions...")
    merged_regions = merge_blast_hits(
        filtered_hits,
        merge_dir,
        DEFAULT_BEDTOOLS_PARAMS,
        deps['bedtools'],
        genome_id=genome_id
    )
    
    merged_output = merge_dir / "merged_regions.tsv"
    save_merged_regions(merged_regions, merged_output)
    
    # Score density calculation (silent)
    region_annotations = annotate_regions_with_best_proteins(
        merged_regions,
        filtered_hits
    )
    
    annotations_output = final_dir / "region_annotations.tsv"
    save_region_annotations(region_annotations, annotations_output)
    
    # STEP 3: Pseudogene detection
    logger.info(f"[3/3] Analyzing {len(region_annotations)} regions for pseudogenes...")
    
    pseudogene_annotations = annotate_pseudogenes(
        region_annotations,
        genome_fasta,
        min_disablements
    )
    
    pseudogenes_output = final_dir / "pseudogene_annotations.tsv"
    save_pseudogene_annotations(pseudogene_annotations, pseudogenes_output)
    
    # Generate summary report
    summary = generate_summary_report(pseudogene_annotations)
    print(summary)
    
    report_file = final_dir / "summary_report.txt"
    with open(report_file, 'w') as f:
        f.write(summary)
    
    num_pseudogenes = sum(1 for ann in pseudogene_annotations if ann.is_pseudogene)
    logger.info(f"Complete: {len(pseudogene_annotations)} regions, {num_pseudogenes} pseudogenes")


def main():
    """Função principal do pipeline."""
    # Parse dos argumentos
    args = parse_arguments()
    
    # Configura logging
    log_file = args.log_file if args.log_file else args.output / "pipeline.log"
    setup_logging(log_file, args.verbose)
    
    logger = logging.getLogger(__name__)
    
    try:
        # Verifica dependências
        deps = verify_and_install_dependencies()
        
        # Decide modo de execução: single genome vs. batch
        if args.input_dir:
            # ========== MODO LOTE ==========
            logger.info("Mode: BATCH (multiple genomes)")
            
            from src.batch_processor import BatchProcessor
            
            processor = BatchProcessor(
                input_dir=args.input_dir,
                lysozymes_path=args.lysozymes,
                output_dir=args.output,
                dependencies=deps,
                min_identity=args.min_identity,
                min_score=args.min_score,
                min_disablements=args.min_disablements,
                num_threads=args.num_threads
            )
            
            start_time = datetime.now()
            processor.run_batch()
            end_time = datetime.now()
            
        else:
            # ========== MODO SINGLE GENOME ==========
            logger.info("Mode: SINGLE GENOME")
            
            # Valida entradas
            validate_inputs(args.genome, args.lysozymes)
            
            # Extrai genome_id do nome do arquivo
            genome_id = args.genome.stem
            
            # Executa pipeline
            start_time = datetime.now()
            
            run_pipeline(
                genome_fasta=args.genome,
                lysozyme_fasta=args.lysozymes,
                output_dir=args.output,
                genome_id=genome_id,
                min_identity=args.min_identity,
                min_score=args.min_score,
                min_disablements=args.min_disablements,
                skip_ssearch=args.skip_ssearch,
                num_threads=args.num_threads,
                deps=deps  # Pass deps to avoid re-verification
            )
            
            end_time = datetime.now()
        
        elapsed = end_time - start_time
        
        logger.debug(f"Total execution time: {elapsed}")
        logger.debug("Pipeline completed successfully!")
        
    except Exception as e:
        logger.error(f"Erro durante execução do pipeline: {e}", exc_info=True)
        sys.exit(1)
        
    except Exception as e:
        logger.error(f"Erro durante execução do pipeline: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
