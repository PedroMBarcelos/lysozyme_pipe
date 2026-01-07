#!/usr/bin/env python3
import pandas as pd

# Load and filter
df = pd.read_csv('BVBRC_genome_ALL.csv')
df_filtered = df[(df['Genome Status'] == 'Complete') & (df['Assembly Accession'].notna())].copy()

print(f"Total: {len(df)}, Complete+Assembly: {len(df_filtered)}")

# Search
results = set()
genome_data = {}

# Lab strains
for marker in ['K-12', 'K12', 'MG1655', 'W3110', 'BW25113', 'DH5', 'DH10B', 'BL21', 'HS', 'Crooks', 'REL606']:
    m = df_filtered[
        df_filtered['Genome Name'].str.contains(marker, na=False, case=False) | 
        df_filtered['Strain'].str.contains(marker, na=False, case=False)
    ]
    for _, row in m.iterrows():
        gid = row['Genome ID']
        results.add(gid)
        genome_data[gid] = row

# Build output
out = []
for gid in sorted(results):
    r = genome_data[gid]
    out.append({
        'Genome ID': gid,
        'Name': r['Genome Name'],
        'Strain': r['Strain'],
        'Assembly': r['Assembly Accession'],
        'Pathovar': r['Pathovar']
    })

df_out = pd.DataFrame(out)
print(f"\nFound {len(df_out)} non-pathogenic genomes:\n")
print(df_out.to_string(index=False))

# Save top 10
df_out.head(10).to_csv('temp_nonpath.csv', index=False)
print(f"\nSaved top 10 to temp_nonpath.csv")
