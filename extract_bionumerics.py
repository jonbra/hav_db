# =============================================================================
# Extract sequences and metadata from BioNumerics export
# =============================================================================

import os
import base64
import csv
import re
from xml.etree import ElementTree as ET
from pathlib import Path

# Set paths
export_dir = Path("export")
entries_dir = export_dir / "_Entries"
experiments_dir = export_dir / "_Experiments"
output_fasta = Path("data/bionumerics_sequences.fasta")
output_metadata = Path("data/bionumerics_metadata.csv")

def decode_b64(text):
    """Decode Base64 string to UTF-8 text"""
    if not text:
        return ""
    try:
        return base64.b64decode(text).decode('utf-8')
    except:
        return text

def get_text_decoded(element):
    """Get text content, decoding Base64 if needed"""
    if element is None:
        return ""
    text = element.text or ""
    if element.get("B64") == "1" and text:
        return decode_b64(text)
    return text

# -----------------------------------------------------------------------------
# Parse all Entry files to get metadata
# -----------------------------------------------------------------------------
print("Parsing metadata from _Entries...")

all_entries = []
entry_files = sorted(entries_dir.glob("*.xml"))

for ef in entry_files:
    print(f"  Reading: {ef.name}")
    tree = ET.parse(ef)
    root = tree.getroot()
    
    for entry in root.findall(".//Entry"):
        # Get sample key (ID)
        key_node = entry.find(".//Key")
        sample_id = get_text_decoded(key_node) if key_node is not None else ""
        
        # Get experiment path (links to sequence file)
        exp_path_node = entry.find(".//SEQ/Path")
        exp_path = get_text_decoded(exp_path_node) if exp_path_node is not None else ""
        
        # Get all fields
        entry_data = {
            "sample_id": sample_id,
            "exp_path": exp_path
        }
        
        for field in entry.findall(".//Field"):
            name_node = field.find("Name")
            value_node = field.find("Value")
            
            field_name = get_text_decoded(name_node) if name_node is not None else ""
            field_value = get_text_decoded(value_node) if value_node is not None else ""
            
            entry_data[field_name] = field_value
        
        all_entries.append(entry_data)

print(f"Found {len(all_entries)} entries")

# -----------------------------------------------------------------------------
# Parse Experiment files to get sequences
# -----------------------------------------------------------------------------
print("\nParsing sequences from _Experiments...")

sequences = {}
missing_files = 0

for i, entry in enumerate(all_entries):
    exp_path = entry.get("exp_path", "")
    sample_id = entry.get("sample_id", "")
    
    if not exp_path or not sample_id:
        continue
    
    # Convert Windows path separator
    exp_file = export_dir / exp_path.replace("\\", "/")
    
    if exp_file.exists():
        try:
            tree = ET.parse(exp_file)
            root = tree.getroot()
            seq_node = root.find(".//SequenceData")
            
            if seq_node is not None and seq_node.text:
                seq_data = seq_node.text.strip().upper()
                if seq_data:
                    sequences[sample_id] = seq_data
        except Exception as e:
            print(f"  Error reading {exp_file}: {e}")
    else:
        missing_files += 1
    
    if (i + 1) % 100 == 0:
        print(f"  Processed {i + 1} of {len(all_entries)}")

print(f"Found {len(sequences)} sequences")
if missing_files > 0:
    print(f"  ({missing_files} experiment files not found)")

# -----------------------------------------------------------------------------
# Create output directory
# -----------------------------------------------------------------------------
output_fasta.parent.mkdir(exist_ok=True)

# -----------------------------------------------------------------------------
# Create output metadata (only for samples with sequences)
# -----------------------------------------------------------------------------
print("\nCreating output files...")

def convert_date(date_str):
    """Convert DD.MM.YYYY to YYYY-MM-DD"""
    if not date_str:
        return ""
    match = re.match(r"(\d{2})\.(\d{2})\.(\d{4})", date_str)
    if match:
        day, month, year = match.groups()
        return f"{year}-{month}-{day}"
    return date_str

# Prepare metadata rows
metadata_rows = []
for entry in all_entries:
    sample_id = entry.get("sample_id", "")
    if sample_id not in sequences:
        continue
    
    row = {
        "sample_id": sample_id,
        "sampling_date": convert_date(entry.get("SAMPLE_DATE", "")),
        "sample_year": entry.get("SAMPLE_YEAR", ""),
        "genotype": entry.get("GENOTYPE", ""),
        "outbreak_variant": entry.get("OUTBREAK_VARIANT", ""),
        "origin": entry.get("ORIGIN", ""),
        "source": entry.get("SOURCE", ""),
        "transmission_route": entry.get("TRANSMISSION_ROUTE", ""),
        "patient_initials": entry.get("PATIENT_INITIALS", ""),
        "group": entry.get("GROUP", ""),
        "sequence_length": len(sequences[sample_id])
    }
    metadata_rows.append(row)

# Write metadata CSV
fieldnames = ["sample_id", "sampling_date", "sample_year", "genotype", 
              "outbreak_variant", "origin", "source", "transmission_route",
              "patient_initials", "group", "sequence_length"]

with open(output_metadata, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(metadata_rows)

print(f"Wrote metadata to: {output_metadata}")

# -----------------------------------------------------------------------------
# Write FASTA file
# -----------------------------------------------------------------------------
with open(output_fasta, "w", encoding="utf-8") as f:
    for sample_id, sequence in sequences.items():
        f.write(f">{sample_id}\n")
        f.write(f"{sequence}\n")

print(f"Wrote {len(sequences)} sequences to: {output_fasta}")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
print("\n" + "=" * 50)
print("SUMMARY")
print("=" * 50)
print(f"Total entries in database: {len(all_entries)}")
print(f"Sequences extracted: {len(sequences)}")
print(f"\nOutput files:")
print(f"  - FASTA: {output_fasta}")
print(f"  - Metadata: {output_metadata}")

# Genotype distribution
print("\nGenotype distribution:")
genotypes = {}
for row in metadata_rows:
    gt = row["genotype"] or "(empty)"
    genotypes[gt] = genotypes.get(gt, 0) + 1
for gt, count in sorted(genotypes.items()):
    print(f"  {gt}: {count}")

# Year distribution
print("\nYear distribution:")
years = {}
for row in metadata_rows:
    yr = row["sample_year"] or "(empty)"
    years[yr] = years.get(yr, 0) + 1
for yr, count in sorted(years.items()):
    print(f"  {yr}: {count}")

# Sequence length stats
lengths = [row["sequence_length"] for row in metadata_rows]
if lengths:
    print(f"\nSequence lengths:")
    print(f"  Min: {min(lengths)} bp")
    print(f"  Max: {max(lengths)} bp")
    print(f"  Mean: {sum(lengths) // len(lengths)} bp")
