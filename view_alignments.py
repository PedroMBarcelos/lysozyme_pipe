import pandas as pd
import sys
import textwrap

def print_alignment(row):
    print("="*80)
    print(f"Region: {row['chromosome']}:{row['start']}-{row['end']}")
    print(f"Protein: {row['best_protein_id']}")
    print(f"Classification: {row.get('classification', 'N/A')}")
    print("-" * 80)
    
    # Get sequences
    qseq = row['qseq'] # Reference Protein
    sseq = row['sseq'] # Genomic Sequence (Translated)
    
    # Handle NaN or non-string values
    if pd.isna(qseq): qseq = ""
    if pd.isna(sseq): sseq = ""
    
    # Calculate match line (identity)
    match_line = ""
    for q, s in zip(qseq, sseq):
        if q == s:
            match_line += "|"
        elif q == '-' or s == '-':
            match_line += " "
        elif s == '*':
            match_line += "!" # Highlight stop codons
        else:
            match_line += "." # Mismatch
            
    # Wrap text for readability (60 chars per line)
    width = 60
    if len(qseq) > 0:
        for i in range(0, len(qseq), width):
            print(f"Ref:  {qseq[i:i+width]}")
            print(f"      {match_line[i:i+width]}")
            print(f"Gen:  {sseq[i:i+width]}")
            print("")
    else:
        print("(No alignment sequence available)")
    print("\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python view_alignments.py <path_to_pseudogene_annotations.tsv>")
        sys.exit(1)
        
    file_path = sys.argv[1]
    
    try:
        # Load the TSV
        df = pd.read_csv(file_path, sep='\t')
        
        # Check if required columns exist
        if 'qseq' not in df.columns or 'sseq' not in df.columns:
            print("Error: The file does not contain 'qseq' and 'sseq' columns.")
            print("Available columns:", df.columns.tolist())
            sys.exit(1)
            
        # Print alignments
        for index, row in df.iterrows():
            print_alignment(row)
            
    except Exception as e:
        print(f"Error reading file: {e}")
