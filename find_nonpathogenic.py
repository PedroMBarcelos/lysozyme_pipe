#!/usr/bin/env python3
"""
Script to find non-pathogenic E. coli genomes from BVBRC database.
Searches for lab strains, commensal isolates, and explicitly marked non-pathogenic genomes.
"""

import pandas as pd
import re

def find_nonpathogenic_genomes(csv_file='BVBRC_genome_ALL.csv'):
    """
    Search for non-pathogenic E. coli genomes using multiple criteria.
    
    Criteria:
    1. Pathovar column explicitly marked as NOT/NON/COMMENSAL
    2. Lab strain markers in Genome Name or Strain fields
    3. Isolation source containing 'laboratory' or 'commensal'
    
    Only includes genomes with:
    - Genome Status = Complete
    - Assembly Accession available (not null)
    """
    
    print("="*80)
    print("SEARCHING FOR NON-PATHOGENIC E. COLI GENOMES IN BVBRC DATABASE")
    print("="*80)
    
    # Load the full database
    df = pd.read_csv(csv_file)
    print(f"\nTotal genomes in database: {len(df)}")
    
    # Filter for Complete genomes with Assembly Accession
    df_filtered = df[
        (df['Genome Status'] == 'Complete') & 
        (df['Assembly Accession'].notna())
    ].copy()
    print(f"Complete genomes with Assembly Accession: {len(df_filtered)}")
    
    # Initialize results set to avoid duplicates
    nonpathogenic_genomes = set()
    genome_data = {}
    
    # Criterion 1: Explicit pathovar marking
    print("\n" + "-"*80)
    print("CRITERION 1: Explicit Pathovar Marking (NOT/NON/COMMENSAL)")
    print("-"*80)
    
    pathovar_patterns = ['NOT', 'NON', 'COMMENSAL', 'NON-PATHOGENIC', 'NONPATHOGENIC']
    for pattern in pathovar_patterns:
        matches = df_filtered[df_filtered['Pathovar'].str.upper().str.contains(pattern, na=False, case=False)]
        if len(matches) > 0:
            print(f"  Pattern '{pattern}': {len(matches)} genomes")
            for _, row in matches.iterrows():
                gid = row['Genome ID']
                nonpathogenic_genomes.add(gid)
                genome_data[gid] = row
    
    # Criterion 2: Lab strain markers
    print("\n" + "-"*80)
    print("CRITERION 2: Laboratory Strain Markers")
    print("-"*80)
    
    lab_markers = [
        'K-12', 'K12', 'MG1655', 'W3110', 'BW25113', 'BW2952', 
        'DH5', 'DH10B', 'BL21', 'HS', 'Crooks', 'REL606',
        'B REL606', 'C ATCC', 'W ATCC'
    ]
    
    for marker in lab_markers:
        # Search in Genome Name
        name_matches = df_filtered[df_filtered['Genome Name'].str.contains(marker, na=False, case=False)]
        # Search in Strain
        strain_matches = df_filtered[df_filtered['Strain'].str.contains(marker, na=False, case=False)]
        
        matches = pd.concat([name_matches, strain_matches]).drop_duplicates(subset=['Genome ID'])
        
        if len(matches) > 0:
            print(f"  Marker '{marker}': {len(matches)} genomes")
            for _, row in matches.iterrows():
                gid = row['Genome ID']
                nonpathogenic_genomes.add(gid)
                genome_data[gid] = row
    
    # Criterion 3: Isolation source
    print("\n" + "-"*80)
    print("CRITERION 3: Isolation Source (laboratory/commensal)")
    print("-"*80)
    
    source_patterns = ['laboratory', 'lab strain', 'commensal', 'non-pathogenic']
    for pattern in source_patterns:
        matches = df_filtered[df_filtered['Isolation Source'].str.contains(pattern, na=False, case=False)]
        if len(matches) > 0:
            print(f"  Pattern '{pattern}': {len(matches)} genomes")
            for _, row in matches.iterrows():
                gid = row['Genome ID']
                nonpathogenic_genomes.add(gid)
                genome_data[gid] = row
    
    # Compile results
    print("\n" + "="*80)
    print(f"TOTAL NON-PATHOGENIC GENOMES FOUND: {len(nonpathogenic_genomes)}")
    print("="*80)
    
    # Create results DataFrame
    results = []
    for gid in sorted(nonpathogenic_genomes):
        row = genome_data[gid]
        results.append({
            'Genome ID': gid,
            'Genome Name': row['Genome Name'],
            'Strain': row['Strain'],
            'Assembly Accession': row['Assembly Accession'],
            'Pathovar': row['Pathovar'],
            'Isolation Source': row['Isolation Source']
        })
    
    results_df = pd.DataFrame(results)
    
    # Display all results
    print("\nFull list of non-pathogenic genomes:\n")
    print(results_df.to_string(index=False))
    
    return results_df

def update_csv_files(selected_genomes_df, output_file='BVBRC_nonpathogenic_genomes.csv'):
    """
    Update the non-pathogenic genomes CSV and the combined CSV.
    """
    
    # Load the full database to get all columns
    full_df = pd.read_csv('BVBRC_genome_ALL.csv')
    
    # Filter full database for selected genome IDs
    selected_ids = selected_genomes_df['Genome ID'].tolist()
    new_nonpath = full_df[full_df['Genome ID'].isin(selected_ids)].copy()
    
    # Ensure Pathovar is set to 'NOT'
    new_nonpath['Pathovar'] = 'NOT'
    
    # Save non-pathogenic genomes file
    new_nonpath.to_csv(output_file, index=False)
    print(f"\n✓ Saved {len(new_nonpath)} non-pathogenic genomes to {output_file}")
    
    # Update combined file
    pathogenic_df = pd.read_csv('BVBRC_genome (1).csv')
    combined_df = pd.concat([pathogenic_df, new_nonpath], ignore_index=True)
    combined_df.to_csv('BVBRC_genome_combined.csv', index=False)
    print(f"✓ Updated BVBRC_genome_combined.csv with {len(combined_df)} total genomes")
    print(f"  - Pathogenic: {len(pathogenic_df)}")
    print(f"  - Non-pathogenic: {len(new_nonpath)}")
    
    return combined_df

if __name__ == '__main__':
    # Find all non-pathogenic genomes
    results_df = find_nonpathogenic_genomes()
    
    if len(results_df) == 0:
        print("\n⚠ No non-pathogenic genomes found!")
        exit(1)
    
    # Ask user how many to add
    print("\n" + "="*80)
    print("SELECTION")
    print("="*80)
    print(f"\nFound {len(results_df)} non-pathogenic genomes.")
    print("\nOptions:")
    print("  1. Add all genomes")
    print("  2. Add specific number (e.g., 10, 20)")
    print("  3. Cancel")
    
    choice = input("\nYour choice (1/2/3): ").strip()
    
    if choice == '1':
        selected_df = results_df
        print(f"\n✓ Selected all {len(selected_df)} genomes")
    elif choice == '2':
        num = int(input("How many genomes to add? "))
        if num > len(results_df):
            print(f"⚠ Requested {num} but only {len(results_df)} available. Adding all.")
            selected_df = results_df
        else:
            selected_df = results_df.head(num)
            print(f"\n✓ Selected first {len(selected_df)} genomes")
    else:
        print("\n✗ Cancelled")
        exit(0)
    
    # Update CSV files
    print("\n" + "="*80)
    print("UPDATING CSV FILES")
    print("="*80)
    
    combined_df = update_csv_files(selected_df)
    
    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
    print("\nNext steps:")
    print("  1. Review the updated CSV files")
    print("  2. Download genomes: python download_genomes.py BVBRC_genome_combined.csv downloaded_genomes/")
    print("  3. Validate downloads: python validate_genomes.py downloaded_genomes/fasta/")
    print("  4. Run pipeline: python pipeline.py --input-dir downloaded_genomes/fasta/ --metadata downloaded_genomes/metadata.tsv")
