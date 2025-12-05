"""
Módulo para detecção de mutações e caracterização de pseudogenes.
Implementa a etapa 6 do pipeline: Detecção de Mutação e Pseudogenes.
"""

import logging
from typing import List, Set, Dict, Optional
from dataclasses import dataclass
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

from src.config import (
    DisablementCounts,
    GAP_CHAR,
    STOP_CODON_CHAR,
    START_CODON_AA,
    MIN_LENGTH_RATIO,
    MAX_LENGTH_RATIO
)
from src.blast_filter import BlastHit
from src.score_density import RegionAnnotation, ProteinHitGroup


logger = logging.getLogger(__name__)


@dataclass
class PseudogeneAnnotation:
    """
    Anotação completa de um pseudogene potencial.
    
    Estrutura de acesso aos atributos:
    - annotation.is_pseudogene
    - annotation.region_annotation.region.chromosome
    - annotation.region_annotation.region.start
    - annotation.region_annotation.region.end
    - annotation.region_annotation.region.strand
    - annotation.region_annotation.best_protein.protein_id
    - annotation.region_annotation.best_protein.score_density
    - annotation.disablements.non_synonymous_substitutions
    - annotation.disablements.in_frame_indels
    - annotation.disablements.frameshifts
    - annotation.disablements.premature_stop_codons
    - annotation.disablements.missing_start_codon
    - annotation.disablements.missing_stop_codon
    - annotation.disablements.total_disablements (propriedade)
    """
    
    region_annotation: RegionAnnotation  # Anotação da região
    disablements: DisablementCounts      # Contadores de mutações
    is_pseudogene: bool                  # Se é classificado como pseudogene
    
    def to_dict(self):
        """Converte a anotação em dicionário."""
        result = self.region_annotation.to_dict()
        result.update(self.disablements.to_dict())
        result['is_pseudogene'] = self.is_pseudogene
        return result


def count_non_synonymous_substitutions(query_seq: str, subject_seq: str) -> int:
    """
    Conta substituições não-sinônimas entre sequências alinhadas.
    
    Conta aminoácidos não idênticos na mesma posição do alinhamento,
    excluindo gaps. IMPORTANTE: Nem toda substituição indica pseudogene,
    apenas contamos para análise.
    
    Args:
        query_seq: Sequência query alinhada (proteína funcional)
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        Número de substituições não-sinônimas
    """
    if len(query_seq) != len(subject_seq):
        logger.warning(
            f"Sequências com comprimentos diferentes: "
            f"query={len(query_seq)}, subject={len(subject_seq)}"
        )
    
    count = 0
    min_length = min(len(query_seq), len(subject_seq))
    
    for i in range(min_length):
        q_aa = query_seq[i]
        s_aa = subject_seq[i]
        
        # Ignora posições com gaps
        if q_aa == GAP_CHAR or s_aa == GAP_CHAR:
            continue
        
        # Ignora stop codons na contagem de substituições
        if q_aa == STOP_CODON_CHAR or s_aa == STOP_CODON_CHAR:
            continue
        
        # Conta se forem diferentes
        if q_aa != s_aa:
            count += 1
    
    # Calcula percentual de diferença
    total_positions = min_length - query_seq.count(GAP_CHAR) - subject_seq.count(GAP_CHAR)
    if total_positions > 0:
        percent_diff = (count / total_positions) * 100
        logger.debug(f"Substituições: {count}/{total_positions} ({percent_diff:.1f}%)")
    
    return count


def count_in_frame_indels(query_seq: str, subject_seq: str) -> int:
    """
    Conta indels in-frame (gaps individuais) presentes em ambas as sequências.
    
    Args:
        query_seq: Sequência query alinhada
        subject_seq: Sequência subject alinhada
    
    Returns:
        Número total de gaps individuais
    """
    query_gaps = query_seq.count(GAP_CHAR)
    subject_gaps = subject_seq.count(GAP_CHAR)
    
    total_gaps = query_gaps + subject_gaps
    
    return total_gaps


def count_frameshifts(hsps: List[BlastHit]) -> int:
    """
    Calcula o número de frameshifts baseado em frames distintos dos HSPs.
    
    Frameshifts = count(distinct_frames) - 1
    
    Args:
        hsps: Lista de HSPs da mesma sequência de referência
    
    Returns:
        Número de frameshifts
    """
    if not hsps:
        return 0
    
    # Coleta todos os frames únicos
    frames: Set[int] = set()
    for hsp in hsps:
        frames.add(hsp.sframe)
    
    # Número de frameshifts é o número de frames distintos menos 1
    num_frameshifts = len(frames) - 1
    
    logger.debug(
        f"Frames detectados: {sorted(frames)}, "
        f"Frameshifts: {num_frameshifts}"
    )
    
    return max(0, num_frameshifts)


def check_missing_start_codon(subject_seq: str) -> bool:
    """
    Verifica se o start codon (Metionina) está ausente.
    
    IMPORTANTE: tblastn traduz o DNA mas NÃO inclui o start codon (ATG->M) 
    automaticamente na sequência traduzida. Por isso, esta verificação é
    relevante apenas se o alinhamento começa no início real do gene.
    
    Args:
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        True se o start codon está ausente, False caso contrário
    """
    # Remove gaps do início
    trimmed_seq = subject_seq.lstrip(GAP_CHAR)
    
    if not trimmed_seq:
        return True
    
    # Para tblastn, a ausência de M no início pode ser normal se o alinhamento
    # não começa exatamente no start codon. Portanto, NÃO consideramos isso
    # como evidência de pseudogene por padrão.
    # Apenas reportamos se estiver completamente ausente E tiver outros indicadores
    first_aa = trimmed_seq[0]
    
    # Relaxa: só marca como problemático se não for M E a sequência for longa o suficiente
    # para esperar que comece no início real
    if len(trimmed_seq) < 50:  # Alinhamentos curtos podem não incluir o start
        return False
    
    is_missing = first_aa != START_CODON_AA
    
    if is_missing:
        logger.debug(f"Possível start codon ausente. Primeiro AA: {first_aa}")
    
    # NÃO marca como problema por padrão - retorna False
    return False


def check_missing_stop_codon(subject_seq: str) -> bool:
    """
    Verifica se o stop codon está ausente no final da sequência.
    
    IMPORTANTE: tblastn traduz até encontrar um stop codon, mas NÃO inclui
    o stop codon (*) na sequência traduzida resultante. Portanto, a ausência
    de * no final é NORMAL e NÃO indica pseudogene.
    
    Args:
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        True se o stop codon está ausente, False caso contrário
    """
    # Remove gaps do final
    trimmed_seq = subject_seq.rstrip(GAP_CHAR)
    
    if not trimmed_seq:
        return True
    
    # Para tblastn, o stop codon NÃO aparece na sequência traduzida normalmente.
    # A presença de * no final seria anormal. A ausência é esperada.
    # Portanto, SEMPRE retornamos False (não é problema)
    last_aa = trimmed_seq[-1]
    
    # Se houver * no final, isso seria estranho, mas não necessariamente problema
    if last_aa == STOP_CODON_CHAR:
        logger.debug(f"Stop codon presente no final (incomum mas não problemático): {last_aa}")
    
    # NÃO marca como problema - retorna sempre False
    return False


def count_premature_stop_codons(subject_seq: str) -> int:
    """
    Conta stop codons prematuros (internos) na sequência.
    
    Args:
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        Número de stop codons internos
    """
    # Remove gaps
    clean_seq = subject_seq.replace(GAP_CHAR, '')
    
    if not clean_seq:
        return 0
    
    # Remove o último caractere (que pode ser stop codon legítimo)
    internal_seq = clean_seq[:-1] if len(clean_seq) > 1 else ""
    
    # Conta stop codons internos
    count = internal_seq.count(STOP_CODON_CHAR)
    
    if count > 0:
        logger.debug(f"Stop codons prematuros detectados: {count}")
    
    return count


def check_length_viability(
    candidate_length: int,
    reference_length: int
) -> Dict[str, any]:
    """
    Apply the 20% Rule to validate gene functionality based on size.
    
    Validates whether candidate gene length is within acceptable range
    (80-120%) of reference protein length. Genes outside this range are
    likely truncated, elongated, or fusion artifacts.
    
    Args:
        candidate_length: Length of candidate sequence in amino acids
        reference_length: Length of reference protein in amino acids
    
    Returns:
        Dict with validation results:
            - candidate_len: Candidate length
            - ref_len: Reference length
            - ratio: Length ratio (candidate/reference)
            - status: "PASS" or "FAIL"
            - classification: Classification result
    """
    # Calculate thresholds (20% Rule)
    min_len = reference_length * MIN_LENGTH_RATIO
    max_len = reference_length * MAX_LENGTH_RATIO
    
    ratio = candidate_length / reference_length if reference_length > 0 else 0
    
    result = {
        "candidate_len": candidate_length,
        "ref_len": reference_length,
        "ratio": round(ratio, 2)
    }
    
    if min_len <= candidate_length <= max_len:
        result["status"] = "PASS"
        result["classification"] = "Size compatible"
    elif candidate_length < min_len:
        result["status"] = "FAIL"
        result["classification"] = "Truncated (size < 80%)"
    else:
        result["status"] = "FAIL"
        result["classification"] = "Elongated/Fusion (size > 120%)"
    
    return result


def check_start_stop_codons_genomic(
    genome_fasta_path: Path,
    chromosome: str,
    start: int,
    end: int,
    strand: str,
    padding: int = 150,
    check_start: bool = True,
    check_stop: bool = True
) -> Dict[str, bool]:
    """
    Check for start/stop codons by expanding aligned region in genome.
    
    Expand & Check strategy with flexible frame verification:
    1. Extract genomic region with padding
    2. Orient sequence (reverse complement if negative strand)
    3. Search for start codon upstream (if check_start=True)
    4. Search for stop codon downstream (if check_stop=True)
    
    Args:
        genome_fasta_path: Path to genome FASTA file
        chromosome: Chromosome/contig ID
        start: Alignment start (1-based)
        end: Alignment end (1-based)
        strand: '+' or '-'
        padding: Extra nucleotides to search (default: 150bp = 50 codons)
        check_start: Search for start codon (default: True)
        check_stop: Search for stop codon (default: True)
    
    Returns:
        Dict with 'has_start', 'has_stop', 'start_codon', 'stop_codon'
    """
    # Load genome sequence
    genome_seqs = {rec.id: rec.seq for rec in SeqIO.parse(genome_fasta_path, "fasta")}
    
    if chromosome not in genome_seqs:
        logger.warning(f"Chromosome {chromosome} not found in genome")
        return {"has_start": False, "has_stop": False, "start_codon": None, "stop_codon": None}
    
    genome_seq = genome_seqs[chromosome]
    
    # Convert to 0-based
    start_0 = start - 1
    end_0 = end
    
    # Adjust coordinates with padding
    search_start = max(0, start_0 - padding)
    search_end = min(len(genome_seq), end_0 + padding)
    
    # Extract extended region
    raw_seq = genome_seq[search_start:search_end]
    
    # Orient sequence (reverse complement if negative strand)
    if strand == '-':
        raw_seq = raw_seq.reverse_complement()
    
    dna_str = str(raw_seq).upper()
    
    # Calculate relative alignment position within extracted string
    if strand == '+':
        rel_align_start = start_0 - search_start
        rel_align_end = end_0 - search_start
    else:
        # On reverse strand, after reverse_complement: original 'end' becomes 5' start
        rel_align_start = search_end - end_0
        rel_align_end = search_end - start_0
    
    # Valid bacterial start codons (E. coli)
    valid_starts = ["ATG", "GTG", "TTG"]
    stop_codons = ["TAA", "TAG", "TGA"]
    
    # --- START CODON SEARCH (UPSTREAM) ---
    has_start = False
    start_codon = None
    
    if check_start:
        upstream_region = dna_str[0:rel_align_start]
        
        # Iterate backwards searching for any start codon (no frame restriction)
        for i in range(len(upstream_region) - 3, -1, -1):
            codon = upstream_region[i:i+3]
            if len(codon) < 3:
                continue
            
            if codon in valid_starts:
                has_start = True
                start_codon = codon
                break
            
            if codon in stop_codons:
                # Stop before start = truncated gene
                break
    
    # --- STOP CODON SEARCH (DOWNSTREAM) ---
    has_stop = False
    stop_codon = None
    
    if check_stop:
        downstream_region = dna_str[rel_align_end:]
        
        # Iterate forward searching for any stop codon (no frame restriction)
        for i in range(0, len(downstream_region) - 2):
            codon = downstream_region[i:i+3]
            if len(codon) < 3:
                continue
            
            if codon in stop_codons:
                has_stop = True
                stop_codon = codon
                break
    
    return {
        "has_start": has_start,
        "has_stop": has_stop,
        "start_codon": start_codon,
        "stop_codon": stop_codon
    }


def analyze_protein_for_disablements(
    protein_hit_group: ProteinHitGroup
) -> DisablementCounts:
    """
    Analyze HSP group for inactivating mutations.
    
    Args:
        protein_hit_group: Protein HSP group
    
    Returns:
        Counters for different mutation types
    """
    logger.debug(f"Analyzing protein {protein_hit_group.protein_id}")
    
    counts = DisablementCounts()
    
    hsps = protein_hit_group.hsps
    
    if not hsps:
        return counts
    
    # Concatenate all HSP sequences for global analysis
    all_query_seq = ""
    all_subject_seq = ""
    
    for hsp in hsps:
        all_query_seq += hsp.qseq
        all_subject_seq += hsp.sseq
    
    # 1. Non-synonymous substitutions
    counts.non_synonymous_substitutions = count_non_synonymous_substitutions(
        all_query_seq, all_subject_seq
    )
    
    # 2. In-frame indels
    counts.in_frame_indels = count_in_frame_indels(
        all_query_seq, all_subject_seq
    )
    
    # 3. Frameshifts
    counts.frameshifts = count_frameshifts(hsps)
    
    # 4. Missing start codon (will be verified genomically in annotate_pseudogenes)
    counts.missing_start_codon = 1 if check_missing_start_codon(all_subject_seq) else 0
    
    # 5. Missing stop codon (will be verified genomically in annotate_pseudogenes)
    counts.missing_stop_codon = 1 if check_missing_stop_codon(all_subject_seq) else 0
    
    # 6. Premature stop codons
    counts.premature_stop_codons = count_premature_stop_codons(all_subject_seq)
    
    logger.debug(f"{protein_hit_group.protein_id}: {counts.total_disablements} mutations")
    
    return counts


def classify_as_pseudogene(
    disablements: DisablementCounts,
    min_disablements: int = 1
) -> bool:
    """
    Classifica se uma região é um pseudogene baseado no número de mutações.
    
    CRITÉRIOS ATUALIZADOS:
    - Stop codons prematuros são forte evidência de pseudogene
    - Frameshifts são forte evidência de pseudogene
    - Ausência de start codon (verificada genomicamente) é forte evidência
    - Ausência de stop codon (verificada genomicamente) é forte evidência
    - Tamanho <80% ou >120% da proteína de referência é forte evidência (Regra dos 20%)
    - Substituições e indels sozinhos NÃO indicam pseudogene (podem ser variação normal)
    
    Args:
        disablements: Contadores de mutações
        min_disablements: Número mínimo de mutações para classificar como pseudogene
    
    Returns:
        True se for classificado como pseudogene, False caso contrário
    """
    # Evidências fortes de pseudogene
    strong_evidence = (
        disablements.premature_stop_codons +
        disablements.frameshifts +
        disablements.missing_start_codon +
        disablements.missing_stop_codon +
        disablements.size_mismatch
    )
    
    # Classifica como pseudogene apenas se houver evidência forte
    is_pseudogene = strong_evidence >= min_disablements
    
    logger.debug(
        f"Classificação: strong_evidence={strong_evidence}, "
        f"is_pseudogene={is_pseudogene}"
    )
    
    return is_pseudogene


def annotate_pseudogenes(
    region_annotations: List[RegionAnnotation],
    genome_fasta_path: Path,
    min_disablements: int = 1,
    padding: int = 150
) -> List[PseudogeneAnnotation]:
    """
    Anota regiões genômicas como potenciais pseudogenes.
    
    Usa verificação genômica "Expand & Check" para detectar ausência
    de start/stop codons na sequência de DNA real.
    
    Args:
        region_annotations: Lista de anotações de região com melhores proteínas
        genome_fasta_path: Caminho para o arquivo FASTA do genoma
        min_disablements: Número mínimo de mutações para classificar como pseudogene
        padding: Nucleotídeos extras para buscar start/stop (padrão: 150bp)
    
    Returns:
        Lista de anotações de pseudogenes
    """
    logger.debug(f"Detecting pseudogenes in {len(region_annotations)} regions...")
    
    pseudogene_annotations = []
    
    for region_ann in region_annotations:
        # Analyze best protein for inactivating mutations
        disablements = analyze_protein_for_disablements(region_ann.best_protein)
        
        # --- SMART START/STOP VERIFICATION ---
        # Step 1: Quick check in alignment (if Methionine present, start codon exists)
        # Step 2: Genomic verification only if alignment check inconclusive
        
        hsps = region_ann.best_protein.hsps
        sorted_hsps = sorted(hsps, key=lambda h: h.qstart)
        first_hsp = sorted_hsps[0]
        last_hsp = sorted_hsps[-1]
        
        # Check start: if alignment has Methionine at beginning (any qstart), start codon exists
        first_qseq_clean = first_hsp.qseq.lstrip('-')
        has_start_in_alignment = (
            len(first_qseq_clean) > 0 and 
            first_qseq_clean[0] == 'M'
            # Note: qstart position doesn't matter - divergent N-terminals are normal
        )
        
        # Check stop: if alignment reaches end of query, stop codon exists
        has_stop_in_alignment = (last_hsp.qend == last_hsp.qlen)
        
        # Genomic verification (only if alignment check failed)
        need_genomic_check = not has_start_in_alignment or not has_stop_in_alignment
        
        if need_genomic_check:
            genomic_check = check_start_stop_codons_genomic(
                genome_fasta_path,
                region_ann.region.chromosome,
                region_ann.region.start,
                region_ann.region.end,
                region_ann.region.strand,
                padding=padding,
                check_start=not has_start_in_alignment,
                check_stop=not has_stop_in_alignment
            )
            
            # Update based on genomic check (only if alignment didn't confirm)
            if not has_start_in_alignment:
                disablements.missing_start_codon = 0 if genomic_check["has_start"] else 1
            else:
                disablements.missing_start_codon = 0  # Confirmed in alignment
            
            if not has_stop_in_alignment:
                disablements.missing_stop_codon = 0 if genomic_check["has_stop"] else 1
            else:
                disablements.missing_stop_codon = 0  # Confirmed in alignment
        else:
            # Both confirmed in alignment - no genomic search needed
            disablements.missing_start_codon = 0
            disablements.missing_stop_codon = 0
        
        # --- SIZE-BASED PSEUDOGENE DETECTION (20% RULE) ---
        # TEMPORARILY DISABLED - Under investigation
        # TODO: Review biological assumptions and implementation
        """
        # Apply only if both start and stop codons are present (per specification)
        has_both_codons = (disablements.missing_start_codon == 0 and 
                          disablements.missing_stop_codon == 0)
        
        if has_both_codons:
            # Calculate candidate length from query coverage (handles overlapping HSPs)
            hsps = region_ann.best_protein.hsps
            min_qstart = min(hsp.qstart for hsp in hsps)
            max_qend = max(hsp.qend for hsp in hsps)
            candidate_length_aa = max_qend - min_qstart + 1
            
            # Reference protein length (all HSPs share same qlen)
            reference_length_aa = hsps[0].qlen
            length_check = check_length_viability(candidate_length_aa, reference_length_aa)
            
            if length_check["status"] == "FAIL":
                disablements.size_mismatch = 1
                logger.debug(
                    f"Size mismatch: candidate={length_check['candidate_len']}aa, "
                    f"reference={length_check['ref_len']}aa, "
                    f"ratio={length_check['ratio']:.2f} - {length_check['classification']}"
                )
            else:
                disablements.size_mismatch = 0
                logger.debug(
                    f"Size check passed: candidate={length_check['candidate_len']}aa, "
                    f"reference={length_check['ref_len']}aa, ratio={length_check['ratio']:.2f}"
                )
        else:
            # Skip size check if start/stop codons are missing
            disablements.size_mismatch = 0
        """
        # Force size_mismatch to 0 while disabled
        disablements.size_mismatch = 0
        
        # Classify as pseudogene
        is_pseudogene = classify_as_pseudogene(disablements, min_disablements)
        
        pseudogene_ann = PseudogeneAnnotation(
            region_annotation=region_ann,
            disablements=disablements,
            is_pseudogene=is_pseudogene
        )
        
        pseudogene_annotations.append(pseudogene_ann)
    
    num_pseudogenes = sum(1 for ann in pseudogene_annotations if ann.is_pseudogene)
    logger.debug(f"Pseudogenes: {num_pseudogenes}/{len(pseudogene_annotations)}")
    
    return pseudogene_annotations


def save_pseudogene_annotations(
    annotations: List[PseudogeneAnnotation],
    output_path
) -> None:
    """
    Salva anotações de pseudogenes em arquivo TSV.
    
    Args:
        annotations: Lista de anotações de pseudogenes
        output_path: Caminho para o arquivo de saída
    """
    import pandas as pd
    
    logger.debug(f"Saving {len(annotations)} pseudogene annotations to: {output_path}")
    
    data = [ann.to_dict() for ann in annotations]
    df = pd.DataFrame(data)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    
    logger.debug("Pseudogene annotations saved successfully")


def save_coverage_statistics(
    annotations: List[PseudogeneAnnotation],
    output_path
) -> None:
    """
    Save detailed coverage statistics for further analysis.
    
    Creates TSV with coverage ratio, alignment details, and classification
    for discussion about size-based validation approaches.
    
    Args:
        annotations: List of pseudogene annotations
        output_path: Path to coverage statistics output file
    """
    import pandas as pd
    
    logger.debug(f"Saving coverage statistics to: {output_path}")
    
    coverage_data = []
    for ann in annotations:
        hsps = ann.region_annotation.best_protein.hsps
        if hsps:
            min_qstart = min(hsp.qstart for hsp in hsps)
            max_qend = max(hsp.qend for hsp in hsps)
            coverage_len = max_qend - min_qstart + 1
            ref_len = hsps[0].qlen
            coverage_ratio = coverage_len / ref_len if ref_len > 0 else 0
            
            # Calculate genomic length
            genomic_len_nt = ann.region_annotation.region.length
            genomic_len_aa = genomic_len_nt / 3
            
            coverage_data.append({
                'chromosome': ann.region_annotation.region.chromosome,
                'start': ann.region_annotation.region.start,
                'end': ann.region_annotation.region.end,
                'strand': ann.region_annotation.region.strand,
                'protein_id': ann.region_annotation.best_protein.protein_id,
                'is_pseudogene': ann.is_pseudogene,
                'classification': 'Pseudogene' if ann.is_pseudogene else 'Functional',
                'alignment_coverage_aa': coverage_len,
                'reference_length_aa': ref_len,
                'coverage_ratio': round(coverage_ratio, 3),
                'genomic_region_nt': genomic_len_nt,
                'genomic_region_aa': round(genomic_len_aa, 1),
                'num_hsps': len(hsps),
                'qstart_min': min_qstart,
                'qend_max': max_qend,
                'missing_start_codon': ann.disablements.missing_start_codon,
                'missing_stop_codon': ann.disablements.missing_stop_codon,
                'premature_stops': ann.disablements.premature_stop_codons,
                'frameshifts': ann.disablements.frameshifts
            })
    
    df = pd.DataFrame(coverage_data)
    
    # Sort by coverage ratio (ascending) to highlight problematic cases
    df = df.sort_values('coverage_ratio', ascending=True)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    
    logger.debug(f"Coverage statistics saved: {len(coverage_data)} regions")


def generate_summary_report(annotations: List[PseudogeneAnnotation]) -> str:
    """
    Gera um relatório resumido das anotações de pseudogenes.
    
    Args:
        annotations: Lista de anotações de pseudogenes
    
    Returns:
        String com o relatório formatado
    """
    total_regions = len(annotations)
    num_pseudogenes = sum(1 for ann in annotations if ann.is_pseudogene)
    num_functional = total_regions - num_pseudogenes
    
    # Estatísticas de mutações
    total_substitutions = sum(ann.disablements.non_synonymous_substitutions for ann in annotations)
    total_indels = sum(ann.disablements.in_frame_indels for ann in annotations)
    total_frameshifts = sum(ann.disablements.frameshifts for ann in annotations)
    total_missing_start = sum(ann.disablements.missing_start_codon for ann in annotations)
    total_missing_stop = sum(ann.disablements.missing_stop_codon for ann in annotations)
    total_premature_stops = sum(ann.disablements.premature_stop_codons for ann in annotations)
    
    # NOVA SEÇÃO: Análise de cobertura e tamanho
    coverage_stats = []
    for ann in annotations:
        hsps = ann.region_annotation.best_protein.hsps
        if hsps:
            min_qstart = min(hsp.qstart for hsp in hsps)
            max_qend = max(hsp.qend for hsp in hsps)
            coverage_len = max_qend - min_qstart + 1
            ref_len = hsps[0].qlen
            coverage_ratio = coverage_len / ref_len if ref_len > 0 else 0
            
            coverage_stats.append({
                'coverage_len': coverage_len,
                'ref_len': ref_len,
                'ratio': coverage_ratio,
                'is_functional': not ann.is_pseudogene
            })
    
    # Estatísticas de cobertura para genes funcionais
    functional_coverage = [s for s in coverage_stats if s['is_functional']]
    pseudogene_coverage = [s for s in coverage_stats if not s['is_functional']]
    
    # Calcular distribuição de ratios
    if functional_coverage:
        func_ratios = [s['ratio'] for s in functional_coverage]
        func_avg_ratio = sum(func_ratios) / len(func_ratios)
        func_min_ratio = min(func_ratios)
        func_max_ratio = max(func_ratios)
        func_below_80 = sum(1 for r in func_ratios if r < 0.8)
    else:
        func_avg_ratio = func_min_ratio = func_max_ratio = func_below_80 = 0
    
    if pseudogene_coverage:
        pseudo_ratios = [s['ratio'] for s in pseudogene_coverage]
        pseudo_avg_ratio = sum(pseudo_ratios) / len(pseudo_ratios)
    else:
        pseudo_avg_ratio = 0
    
    report = f"""
╭──────────────────────────────────────────────────────────────────╮
│          RELATÓRIO DE ANOTAÇÃO DE PSEUDOGENES DE LISOZIMAS       │
╰──────────────────────────────────────────────────────────────────╯

RESUMO GERAL:
  Total de regiões analisadas:     {total_regions}
  Pseudogenes detectados:          {num_pseudogenes} ({100*num_pseudogenes/total_regions if total_regions > 0 else 0:.1f}%)
  Genes funcionais:                {num_functional} ({100*num_functional/total_regions if total_regions > 0 else 0:.1f}%)

ESTATÍSTICAS DE MUTAÇÕES:
  Substituições não-sinônimas:     {total_substitutions}
  Indels in-frame:                 {total_indels}
  Frameshifts:                     {total_frameshifts}
  Ausência de start codon:         {total_missing_start}
  Ausência de stop codon:          {total_missing_stop}
  Stop codons prematuros:          {total_premature_stops}
  
  Total de mutações inativadoras:  {total_substitutions + total_indels + total_frameshifts + total_missing_start + total_missing_stop + total_premature_stops}

ANÁLISE DE COBERTURA DA PROTEÍNA DE REFERÊNCIA:
  (Razão = Cobertura do Alinhamento / Tamanho da Referência)
  
  Genes Funcionais ({len(functional_coverage)} regiões):
    Razão média de cobertura:      {func_avg_ratio:.2f}
    Razão mínima:                   {func_min_ratio:.2f}
    Razão máxima:                   {func_max_ratio:.2f}
    Regiões com cobertura <80%:     {func_below_80} ({100*func_below_80/len(functional_coverage) if functional_coverage else 0:.1f}%)
  
  Pseudogenes ({len(pseudogene_coverage)} regiões):
    Razão média de cobertura:      {pseudo_avg_ratio:.2f}

  ⚠️  NOTA: Cobertura <80% pode indicar:
     - Alinhamentos parciais (domínio conservado vs proteína completa)
     - Proteínas de referência multi-domínio
     - Genes truncados reais
     ➜ Requer análise caso-a-caso para determinar funcionalidade

"""
    
    return report
