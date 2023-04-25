# Load the required libraries
library(sleuth)
library(dplyr)

# Read in the table describing samples and Kallisto output paths
stab <- read.csv("table.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Initialize the sleuth object using the sleuth_prep function
so <- sleuth_prep(stab)

# Perform differential expression analysis comparing the two conditions
so <- sleuth_fit(so, ~condition, 'full')

# Fit the reduced model to compare in the likelihood ratio test
so <- sleuth_fit(so, ~1, 'reduced')

# Perform the likelihood ratio test for differential expression between conditions
so <- sleuth_lrt(so, 'reduced', 'full')

# Extract the test results from the sleuth object
results_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

# Filter most significant results (FDR/qval < 0.05) and sort by pval
significant_results <- dplyr::filter(results_table, pval <= 0.05) %>% dplyr::arrange(pval)

# Print top 10 transcripts
head(significant_results, n = 10)

#write FDR < 0.05 transcripts to file
write.table(significant_results, file="fdr05_results.txt",quote = FALSE,row.names = FALSE)

# Get the IDs of top transcripts
top_transcripts <- head(select(significant_results, target_id, pval, qval), n = 10)

# Extract needed results from Kallisto for plotting
so <- sleuth_prep(stab, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

# Create a ranked gene list using Ensembl gene IDs
rnk <- data.frame(gene_ids = significant_results$target_id,
                  rank = seq_along(significant_results$target_id))

# Save the ranked gene list as a tab-separated file for WebGestalt
write.table(rnk, file = "gene_ids.rnk", sep = "\t", quote = FALSE, row.names = FALSE)

# Use the first target_id from top_transcripts for plotting
if (nrow(top_transcripts) > 0) {
  top_target_id <- top_transcripts$target_id[1]
  topplot <- plot_bootstrap(so, top_target_id, units = "tpm", color_by = "condition")
  
  # Save plot as PNG file
  png(paste0(top_target_id, "_TPM.png"))
  plot(topplot)
  dev.off()
} else {
  cat("No significant transcripts found.\n")
}
