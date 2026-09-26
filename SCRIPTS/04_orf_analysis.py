# 04_orf_analysis.py
# Detects open reading frames (ORFs) in BRCA1 and BRCA2 across all 3 forward reading frames.
# Assumes brca1_seq and brca2_seq are already loaded (see 01_fetch_and_parse.py).

import pandas as pd
import matplotlib.pyplot as plt


def find_orfs(seq, min_length=100):
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []
    seq_len = len(seq)

    for frame in range(3):
        i = frame
        while i < seq_len - 2:
            codon = seq[i:i + 3]
            if codon == "ATG":
                j = i
                while j < seq_len - 2:
                    stop_codon = seq[j:j + 3]
                    if stop_codon in stop_codons:
                        orf_length = j + 3 - i
                        if orf_length >= min_length:
                            orfs.append({
                                "start": i, "end": j + 3,
                                "length": orf_length, "frame": frame + 1,
                            })
                        i = j
                        break
                    j += 3
                else:
                    break
            i += 3
    return orfs


brca1_orfs = find_orfs(brca1_seq)
brca2_orfs = find_orfs(brca2_seq)

print("BRCA1: Number of ORFs found:", len(brca1_orfs))
print("BRCA2: Number of ORFs found:", len(brca2_orfs))

if brca1_orfs:
    longest_brca1 = max(brca1_orfs, key=lambda x: x["length"])
    print("BRCA1 Longest ORF length:", longest_brca1["length"], "bp")

if brca2_orfs:
    longest_brca2 = max(brca2_orfs, key=lambda x: x["length"])
    print("BRCA2 Longest ORF length:", longest_brca2["length"], "bp")

# Save individual ORF tables
brca1_orf_df = pd.DataFrame(brca1_orfs)
brca2_orf_df = pd.DataFrame(brca2_orfs)
brca1_orf_df.to_csv("results/tables/brca1_orf_summary.csv", index=False)
brca2_orf_df.to_csv("results/tables/brca2_orf_summary.csv", index=False)

# ORF length distribution plot
plt.figure(figsize=(10, 5))
brca1_lengths = [orf["length"] for orf in brca1_orfs]
brca2_lengths = [orf["length"] for orf in brca2_orfs]

plt.hist(brca1_lengths, bins=20, alpha=0.6, label="BRCA1", color="steelblue")
plt.hist(brca2_lengths, bins=20, alpha=0.6, label="BRCA2", color="darkorange")
plt.xlabel("ORF Length (bp)")
plt.ylabel("Frequency")
plt.title("ORF Length Distribution: BRCA1 vs BRCA2")
plt.legend()
plt.grid(alpha=0.3)
plt.savefig("results/figures/orf_length_distribution.png", dpi=300, bbox_inches="tight")
plt.show()
