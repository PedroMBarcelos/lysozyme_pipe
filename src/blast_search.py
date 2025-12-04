"""
Module for BLAST database preparation and search execution.
Implements pipeline step 1: Preparation and Homology Search.
"""

import subprocess
import logging
from pathlib import Path
from typing import Optional

from src.config import BlastParameters, DEFAULT_BLAST_PARAMS


logger = logging.getLogger(__name__)


class BlastDatabaseError(Exception):
    """Exception for BLAST database creation errors."""
    pass


class BlastSearchError(Exception):
    """Exception for BLAST search execution errors."""
    pass


def create_blast_database(
    genome_fasta_path: Path,
    output_db_path: Optional[Path] = None,
    dbtype: str = "nucl",
    makeblastdb_path: str = "makeblastdb"
) -> Path:
    """
    Create a BLAST database from a FASTA file.
    
    Args:
        genome_fasta_path: Path to genome FASTA file
        output_db_path: Path to output database (optional)
        dbtype: Database type ('nucl' for nucleotides, 'prot' for proteins)
    
    Returns:
        Path to created database
    
    Raises:
        BlastDatabaseError: If database creation fails
    """
    if not genome_fasta_path.exists():
        raise BlastDatabaseError(f"Arquivo de genoma não encontrado: {genome_fasta_path}")
    
    # Use same path as input file if not specified
    if output_db_path is None:
        output_db_path = genome_fasta_path
    
    logger.debug(f"Creating BLAST database...")
    
    command = [
        makeblastdb_path,
        "-in", str(genome_fasta_path),
        "-dbtype", dbtype,
        "-out", str(output_db_path)
    ]
    
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )
        logger.debug("BLAST database created")
        return output_db_path
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Erro ao criar banco de dados BLAST: {e.stderr}"
        logger.error(error_msg)
        raise BlastDatabaseError(error_msg) from e
    except FileNotFoundError:
        error_msg = "makeblastdb não encontrado. Certifique-se de que o BLAST+ está instalado."
        logger.error(error_msg)
        raise BlastDatabaseError(error_msg)


def run_tblastn_search(
    query_protein_fasta: Path,
    database_path: Path,
    output_file: Path,
    blast_params: BlastParameters = DEFAULT_BLAST_PARAMS,
    tblastn_path: str = "tblastn"
) -> Path:
    """
    Execute tblastn search with strict parameters.
    
    Search protein sequences (query) against nucleotide database,
    translating database in all 6 reading frames.
    
    Args:
        query_protein_fasta: FASTA file with query protein sequences (lysozymes)
        database_path: Path to nucleotide BLAST database
        output_file: Output file for results
        blast_params: BLAST parameters (uses DEFAULT_BLAST_PARAMS if not specified)
    
    Returns:
        Path to output file with results
    
    Raises:
        BlastSearchError: If search execution fails
    """
    if not query_protein_fasta.exists():
        raise BlastSearchError(f"Arquivo query não encontrado: {query_protein_fasta}")
    
    # Create output directory if it doesn't exist
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    logger.debug(f"Running tblastn with {blast_params.num_threads} threads...")
    
    # Build command with parameters
    command = [
        tblastn_path,
        "-query", str(query_protein_fasta),
        "-db", str(database_path),
        "-out", str(output_file)
    ]
    
    # Adiciona os parâmetros customizados
    command.extend(blast_params.to_command_args())
    
    logger.debug(f"BLAST command: {' '.join(command)}")
    
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )
        
        logger.debug("BLAST search complete")
        return output_file
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Erro ao executar tblastn: {e.stderr}"
        logger.error(error_msg)
        raise BlastSearchError(error_msg) from e
    except FileNotFoundError:
        error_msg = "tblastn não encontrado. Certifique-se de que o BLAST+ está instalado."
        logger.error(error_msg)
        raise BlastSearchError(error_msg)


def run_blast_pipeline_step(
    genome_fasta: Path,
    lysozyme_fasta: Path,
    output_dir: Path,
    blast_params: BlastParameters = DEFAULT_BLAST_PARAMS,
    makeblastdb_path: str = "makeblastdb",
    tblastn_path: str = "tblastn"
) -> Path:
    """
    Execute complete BLAST preparation and search step.
    
    Args:
        genome_fasta: E. coli genome FASTA file
        lysozyme_fasta: Reference lysozyme FASTA file
        output_dir: Output directory for files
        blast_params: BLAST parameters
    
    Returns:
        Path to file with BLAST results
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create BLAST database
    db_path = output_dir / "genome_db"
    create_blast_database(genome_fasta, db_path, makeblastdb_path=makeblastdb_path)
    
    # Execute tblastn search
    blast_output = output_dir / "blast_results.tsv"
    run_tblastn_search(lysozyme_fasta, db_path, blast_output, blast_params, tblastn_path)
    
    return blast_output
