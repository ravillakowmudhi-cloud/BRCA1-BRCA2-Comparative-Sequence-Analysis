# 02_sequence_composition.py
# Calculates nucleotide composition (A/T/G/C %, GC%) for BRCA1 and BRCA2.
# Assumes brca1_seq and brca2_seq are already loaded (see 01_fetch_and_parse.py).

import pandas as pd


def nucleotide_composition(seq):
    length = len(seq)
    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")
    gc_percent = (g + c) / length * 100
    at_percent = (a + t) / length * 100
    return {
        "Length (bp)": length,
        "A (%)": round(a / length * 100, 2),
        "T (%)": round(t / length * 100, 2),
        "G (%)": round(g / length * 100, 2),
        "C (%)": round(c / length * 100, 2),
        "GC (%)": round(gc_percent, 2),
        "AT (%)": round(at_percent, 2),
    }


brca1_comp = nucleotide_composition(brca1_seq)
brca2_comp = nucleotide_composition(brca2_seq)

print("BRCA1 Composition:")
for k, v in brca1_comp.items():
    print(f"  {k}: {v}")

print("\nBRCA2 Composition:")
for k, v in brca2_comp.items():
    print(f"  {k}: {v}")

# Save as a table
comparison_df = pd.DataFrame({"BRCA1": brca1_comp, "BRCA2": brca2_comp})
comparison_df.to_csv("results/tables/nucleotide_composition_summary.csv")

# Grouped bar chart: nucleotide composition comparison
import matplotlib.pyplot as plt
import numpy as np

bases = ["A (%)", "T (%)", "G (%)", "C (%)"]
brca1_values = [brca1_comp[b] for b in bases]
brca2_values = [brca2_comp[b] for b in bases]

x = np.arange(len(bases))
width = 0.35

plt.figure(figsize=(8, 5))
plt.bar(x - width / 2, brca1_values, width, label="BRCA1", color="steelblue")
plt.bar(x + width / 2, brca2_values, width, label="BRCA2", color="darkorange")
plt.xlabel("Nucleotide")
plt.ylabel("Percentage (%)")
plt.title("Nucleotide Composition Comparison: BRCA1 vs BRCA2")
plt.xticks(x, ["A", "T", "G", "C"])
plt.legend()
plt.grid(alpha=0.3, axis="y")
plt.savefig("results/figures/nucleotide_composition_barchart.png", dpi=300, bbox_inches="tight")
plt.show()
