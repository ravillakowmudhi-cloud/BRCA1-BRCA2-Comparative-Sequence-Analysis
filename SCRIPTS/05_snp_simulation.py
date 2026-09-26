# 05_snp_simulation.py
# Simulates 20 random point mutations per gene and classifies each as
# synonymous or non-synonymous using Biopython's translation function.
# Assumes brca1_seq and brca2_seq are already loaded (see 01_fetch_and_parse.py).

import random
import pandas as pd
from Bio.Seq import Seq


def simulate_snps(seq, num_mutations=20, seed=42):
    random.seed(seed)
    bases = ["A", "T", "G", "C"]
    max_pos = (len(seq) // 3) * 3
    results = []

    for _ in range(num_mutations):
        pos = random.randrange(0, max_pos)
        codon_start = (pos // 3) * 3
        original_codon = seq[codon_start:codon_start + 3]

        original_base = seq[pos]
        new_base = random.choice([b for b in bases if b != original_base])

        mutated_codon = (
            original_codon[:pos - codon_start]
            + new_base
            + original_codon[pos - codon_start + 1:]
        )

        original_aa = str(Seq(original_codon).translate())
        mutated_aa = str(Seq(mutated_codon).translate())

        mutation_type = "Synonymous" if original_aa == mutated_aa else "Non-synonymous"

        results.append({
            "position": pos,
            "original_codon": original_codon,
            "mutated_codon": mutated_codon,
            "original_aa": original_aa,
            "mutated_aa": mutated_aa,
            "mutation_type": mutation_type,
        })

    return pd.DataFrame(results)


brca1_snp_df = simulate_snps(brca1_seq)
brca2_snp_df = simulate_snps(brca2_seq)

print("BRCA1:\n", brca1_snp_df["mutation_type"].value_counts())
print("\nBRCA2:\n", brca2_snp_df["mutation_type"].value_counts())

brca1_snp_df.to_csv("results/tables/brca1_snp_simulation.csv", index=False)
brca2_snp_df.to_csv("results/tables/brca2_snp_simulation.csv", index=False)
