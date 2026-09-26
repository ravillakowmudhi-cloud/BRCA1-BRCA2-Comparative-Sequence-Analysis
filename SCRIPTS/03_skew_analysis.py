# 03_skew_analysis.py
# Sliding-window GC content, GC skew, and AT skew for BRCA1 and BRCA2.
# Assumes brca1_seq and brca2_seq are already loaded (see 01_fetch_and_parse.py).

import statistics
import matplotlib.pyplot as plt


def sliding_gc(seq, window=100, step=50):
    gc_values = []
    positions = []
    for i in range(0, len(seq) - window + 1, step):
        window_seq = seq[i:i + window]
        gc = (window_seq.count("G") + window_seq.count("C")) / window * 100
        gc_values.append(gc)
        positions.append(i)
    return positions, gc_values


def sliding_skew(seq, window=100, step=50):
    gc_skew, at_skew, positions = [], [], []
    for i in range(0, len(seq) - window + 1, step):
        window_seq = seq[i:i + window]
        g, c = window_seq.count("G"), window_seq.count("C")
        a, t = window_seq.count("A"), window_seq.count("T")
        gc_s = (g - c) / (g + c) if (g + c) != 0 else 0
        at_s = (a - t) / (a + t) if (a + t) != 0 else 0
        gc_skew.append(gc_s)
        at_skew.append(at_s)
        positions.append(i)
    return positions, gc_skew, at_skew


# --- Sliding-window GC content ---
brca1_positions, brca1_gc = sliding_gc(brca1_seq)
brca2_positions, brca2_gc = sliding_gc(brca2_seq)

plt.figure(figsize=(12, 5))
plt.plot(brca1_positions, brca1_gc, label="BRCA1", color="steelblue")
plt.plot(brca2_positions, brca2_gc, label="BRCA2", color="darkorange")
plt.xlabel("Position along sequence (bp)")
plt.ylabel("GC content (%)")
plt.title("Sliding-Window GC Content: BRCA1 vs BRCA2")
plt.legend()
plt.grid(alpha=0.3)
plt.savefig("results/figures/gc_content_windows.png", dpi=300, bbox_inches="tight")
plt.show()

# --- GC / AT skew ---
brca1_pos_skew, brca1_gc_skew, brca1_at_skew = sliding_skew(brca1_seq)
brca2_pos_skew, brca2_gc_skew, brca2_at_skew = sliding_skew(brca2_seq)

print("BRCA1 Mean GC Skew:", round(statistics.mean(brca1_gc_skew), 4))
print("BRCA1 Mean AT Skew:", round(statistics.mean(brca1_at_skew), 4))
print("BRCA2 Mean GC Skew:", round(statistics.mean(brca2_gc_skew), 4))
print("BRCA2 Mean AT Skew:", round(statistics.mean(brca2_at_skew), 4))

plt.figure(figsize=(12, 5))
plt.plot(brca1_pos_skew, brca1_gc_skew, label="BRCA1 GC Skew", color="steelblue")
plt.plot(brca2_pos_skew, brca2_gc_skew, label="BRCA2 GC Skew", color="darkorange")
plt.axhline(0, color="gray", linestyle="--", linewidth=1)
plt.xlabel("Position along sequence (bp)")
plt.ylabel("GC Skew")
plt.title("GC Skew Comparison: BRCA1 vs BRCA2")
plt.legend()
plt.grid(alpha=0.3)
plt.savefig("results/figures/gc_at_skew_plot.png", dpi=300, bbox_inches="tight")
plt.show()

# --- Export window-level GC content for R statistical testing ---
import pandas as pd

gc_export = pd.DataFrame({
    "gene": ["BRCA1"] * len(brca1_gc) + ["BRCA2"] * len(brca2_gc),
    "position": brca1_positions + brca2_positions,
    "gc_content": brca1_gc + brca2_gc,
})
gc_export.to_csv("results/tables/gc_content_for_r.csv", index=False)
