# 01_fetch_and_parse.py
# Fetches human BRCA1 and BRCA2 reference sequences from NCBI and parses them.

# --- Step 1: Fetch FASTA files from NCBI (run in Colab with '!' prefix, or in a terminal without it) ---
# !curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NM_007294.4&rettype=fasta&retmode=text" -o data/fasta/brca1.fasta
# !curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NM_000059.4&rettype=fasta&retmode=text" -o data/fasta/brca2.fasta

from Bio import SeqIO

# --- Step 2: Parse both sequences ---
brca1_record = SeqIO.read("data/fasta/brca1.fasta", "fasta")
brca1_seq = str(brca1_record.seq)

brca2_record = SeqIO.read("data/fasta/brca2.fasta", "fasta")
brca2_seq = str(brca2_record.seq)

print("BRCA1")
print("ID:", brca1_record.id)
print("Description:", brca1_record.description)
print("Length:", len(brca1_seq), "bp")
print()
print("BRCA2")
print("ID:", brca2_record.id)
print("Description:", brca2_record.description)
print("Length:", len(brca2_seq), "bp")
