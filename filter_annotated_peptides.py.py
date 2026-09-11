# =========================================================
# Required packages (install once if needed)
# =========================================================
# pip install pandas biopython openpyxl

import pandas as pd
from Bio import SeqIO
import re
import os

# =========================================================
# File paths
# =========================================================

base_path = r"C:\Users\manju\OneDrive - Manipal Academy of Higher Education\MSc Project Sem\Mapping_P"

peptide_xls = os.path.join(base_path, "E_histolytica_petide20_02_26.xlsx")
protein_fasta = os.path.join(base_path, "AmoebaDB-68_EhistolyticaHM1IMSS_AnnotatedProteins.fasta")

matched_output = os.path.join(base_path, "E_histolytica_matched_peptides.xlsx")
unmatched_output = os.path.join(base_path, "E_histolytica_unmatched_peptides.xlsx")

# =========================================================
# 1. Read peptide Excel file
# =========================================================

peptides_df = pd.read_excel(peptide_xls)

if "Annotated Sequence" not in peptides_df.columns:
    raise ValueError("Column 'Annotated Sequence' not found in Excel file")

# =========================================================
# 2. Clean peptide sequences (same logic as R)
# =========================================================

def clean_peptide(seq):
    seq = str(seq).upper().strip()
    seq = re.sub(r"\(.*?\)", "", seq)   # remove (modifications)
    seq = re.sub(r"\[.*?\]", "", seq)   # remove [modifications]
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
    return seq

peptides_df["Clean_Peptide"] = peptides_df["Annotated Sequence"].apply(clean_peptide)

# Remove empty peptides
peptides_df = peptides_df[peptides_df["Clean_Peptide"] != ""]

# =========================================================
# 3. Read protein FASTA
# =========================================================

protein_sequences = [
    str(record.seq).upper()
    for record in SeqIO.parse(protein_fasta, "fasta")
]

print("Total proteins loaded:", len(protein_sequences))

# =========================================================
# 4. Peptide-to-protein mapping
# =========================================================

def is_peptide_matched(peptide, protein_seqs):
    return any(peptide in protein for protein in protein_seqs)

match_status = peptides_df["Clean_Peptide"].apply(
    lambda x: is_peptide_matched(x, protein_sequences)
)

# =========================================================
# 5. Split matched / unmatched
# =========================================================

matched_peptides = peptides_df[match_status].copy()
unmatched_peptides = peptides_df[~match_status].copy()

# Remove helper column before export
matched_peptides.drop(columns=["Clean_Peptide"], inplace=True)
unmatched_peptides.drop(columns=["Clean_Peptide"], inplace=True)

# =========================================================
# 6. Write output files
# =========================================================

matched_peptides.to_excel(matched_output, index=False)
unmatched_peptides.to_excel(unmatched_output, index=False)

# =========================================================
# 7. Summary
# =========================================================

print("Mapping completed")
print("Matched peptides   :", len(matched_peptides))
print("Unmatched peptides :", len(unmatched_peptides))
print("Total input rows   :", len(matched_peptides) + len(unmatched_peptides))