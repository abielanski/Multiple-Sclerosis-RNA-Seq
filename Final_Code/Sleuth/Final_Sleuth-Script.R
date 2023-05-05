# Install and Load the required libraries

# Devtools - package to help run R script 
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")
library(sleuth)

# Dplyr - package for data manipulation
install.packages("dplyr")
library(dplyr)

# Biomart - a package that provides an interface to BioMart databases (e.g., Ensembl).
install.packages("biomaRt")
library(biomaRt)

# install.packages("ggplot2")
# a powerful package for creating static, interactive, and animated graphics
library(ggplot2)

# Read in the table describing samples and Kallisto output paths   
stab <- read.csv("table.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")

# Create a unique list of conditions from the 'condition' column of your input data
unique_conditions <- unique(stab$condition)

# Initialize the sleuth object using the sleuth_prep function
so <- sleuth_prep(stab)

# Perform differential expression analysis comparing the conditions
so <- sleuth_fit(so, ~condition, 'full')

# Fit the reduced model to compare in the likelihood ratio test
so <- sleuth_fit(so, ~1, 'reduced')

# Perform the likelihood ratio test for differential expression between conditions
so <- sleuth_lrt(so, 'reduced', 'full')

# Extract the test results from the sleuth object
results_table_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

# Initialize an empty list to store Wald test results
wald_test_results <- list()

# Iterate through all pairs of conditions and perform the Wald test
for (i in 1:(length(unique_conditions) - 1)) {
  for (j in (i + 1):length(unique_conditions)) {
    condition1 <- unique_conditions[i]
    condition2 <- unique_conditions[j]
    
    # Filter samples for the current pair of conditions
    stab_filtered <- stab[stab$condition %in% c(condition1, condition2),]
    
    # Re-run sleuth_prep, sleuth_fit, and Wald test for the current pair of conditions
    so_filtered <- sleuth_prep(stab_filtered)
    so_filtered <- sleuth_fit(so_filtered, ~condition, 'full')
    
    # Perform the Wald test and store the results in the wald_test_results list
    test_name <- paste("wt", condition1, "vs", condition2, sep = "_")
    so_filtered <- sleuth_wt(so_filtered, which_beta = paste0("condition", condition2), 'full')
    results_table_wt <- sleuth_results(so_filtered, paste0("condition", condition2), 'wt', show_all = FALSE)
    significant_results_wt <- dplyr::filter(results_table_wt, pval <= 0.05) %>% dplyr::arrange(pval)
    wald_test_results[[test_name]] <- significant_results_wt
  }
}

# Process and save the results for LRT test comparisons
for (test_name in names(list(lrt = results_table_lrt))) {
  results <- list(lrt = results_table_lrt)[[test_name]]
  
  # Write FDR < 0.05 transcripts to file
  output_file <- paste0("fdr05_results_", test_name, ".txt")
  write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Process and save the results for Wald test comparisons
for (test_name in names(wald_test_results)) {
  results <- wald_test_results[[test_name]]
  
  # Write FDR < 0.05 transcripts to file
  output_file <- paste0("fdr05_results_", test_name, ".txt")
  write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # Read the top transcripts from the saved txt file
  txt_file <- output_file
  fdr05_results <- read.table(txt_file, header = TRUE, sep = "\t")
  
  # Get the Ensembl transcript IDs from the txt file
  transcript_ids <- fdr05_results$target_id
  
  # Set up the BioMart object for the organism of interest (e.g., Mus musculus)
  mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  
  # Convert Ensembl transcript IDs to Ensembl gene IDs
  ensembl_transcript_ids <- transcript_ids
  
  t2g <- getBM(filters = "ensembl_transcript_id_version",
               attributes = c("ensembl_transcript_id_version", "ensembl_gene_id",'external_gene_name'),
               values = ensembl_transcript_ids,
               mart = mart)
  
  # Merge the gene IDs with the transcript IDs and p-values from the txt file
  transcript_ids_gene_ids_pvals <- merge(fdr05_results, t2g, by.x = "target_id", by.y = "ensembl_transcript_id_version")
  
  # Calculate the score and create a new data frame with gene_id and score columns
  rnk <- data.frame(gene_id = transcript_ids_gene_ids_pvals$ensembl_gene_id,
                    score = -log10(transcript_ids_gene_ids_pvals$pval))
  
  # Remove any rows with missing gene IDs or scores
  rnk <- na.omit(rnk)
  
  rnk_file <- paste0("gene_ids_with_scores_", test_name, ".rnk")
  write.table(rnk, file = rnk_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # Ensure there are no spaces around the tab delimiter
  rnk_content <- readLines(rnk_file)
  rnk_content <- gsub("\\s+\\t|\\t\\s+", "\t", rnk_content)
  writeLines(rnk_content, rnk_file)
}


# Use the first target_id from top_transcripts for plotting

# Create PCA plot
pca_plot <- plot_pca(so, color_by = "condition", text_labels = FALSE) +
  theme_bw() +
  labs(title = "PCA Plot")

# Save PCA plot as a .png file
ggsave("pca_plot.png", plot = pca_plot, width = 7, height = 5, dpi = 300)

# Create bootstrap plot
so <- sleuth_prep(stab, extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

# Extract the top transcript ID from the significant_results data frame
top_transcript_id <- significant_results_lrt$target_id[1]

# Call plot_bootstrap function
bootstrap_plot <- plot_bootstrap(so, top_transcript_id, units = "tpm", color_by = "condition")

# Save plot as PNG file
png(filename = paste0(top_transcript_id, "_TPM_bootstrap_plot.png"))
print(bootstrap_plot)
dev.off()

# Create a volcano plot using the significant results of the Wald test for each comparison
for (test_name in c("wt_CFA_vs_EAE_Sephin1", "wt_CFA_vs_EAE_Vehicle", "wt_EAE_Sephin1_vs_EAE_Vehicle")) {
  volcano_data <- combined_results[[test_name]] %>%
    dplyr::mutate(neg_log10_pval = -log10(pval)) %>%
    dplyr::select(target_id, b, neg_log10_pval)
  
  volcano_plot <- ggplot(volcano_data, aes(x = b, y = neg_log10_pval)) +
    geom_point(alpha = 0.5, size = 1) +
    theme_bw() +
    labs(title = paste0("Volcano Plot: ", test_name),
         x = "Log2 Fold Change",
         y = "-Log10 P-value")
  
  # Save volcano plot as a .png file
  ggsave(paste0("volcano_plot_", test_name, ".png"), plot = volcano_plot, width = 7, height = 5, dpi = 300)
}
