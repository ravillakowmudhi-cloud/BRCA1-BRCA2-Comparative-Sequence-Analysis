# 06_statistics.R
# Statistical comparison of window-level GC content between BRCA1 and BRCA2.
# Input: gc_content_for_r.csv (exported from 03_skew_analysis.py)

data <- read.csv("gc_content_for_r.csv")
head(data)
str(data)

# --- Welch two-sample t-test ---
brca1_data <- data[data$gene == "BRCA1", "gc_content"]
brca2_data <- data[data$gene == "BRCA2", "gc_content"]

t_test_result <- t.test(brca1_data, brca2_data)
print(t_test_result)

# Save results
capture.output(t_test_result, file = "statistical_test_results.txt")

# --- Boxplot comparison ---
install.packages("ggplot2")
library(ggplot2)

ggplot(data, aes(x = gene, y = gc_content, fill = gene)) +
  geom_boxplot() +
  labs(
    title = "GC Content Comparison: BRCA1 vs BRCA2",
    x = "Gene", y = "GC Content (%)"
  ) +
  theme_minimal()

ggsave("brca1_vs_brca2_boxplot.png", width = 8, height = 5, dpi = 300)
