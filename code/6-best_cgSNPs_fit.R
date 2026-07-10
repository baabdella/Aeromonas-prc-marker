#!/usr/bin/env Rscript

# ==============================================================================
# Phylogenetic Congruence Analysis: Gene Trees vs. cgSNP Reference
# Using ONLY the 'ape' package for Symmetric Difference
# ==============================================================================

if (!requireNamespace("ape", quietly = TRUE)) {
    install.packages("ape", repos = "https://cloud.r-project.org")
}

library(ape)

# --- 1. Configuration ---
ref_tree_path <- "cgSNP_reference.nwk"
gene_tree_folder <- "gene_trees_nwk"
output_csv <- "congruence_results_final.csv"

# --- 2. Helper Functions ---

# Sanitizes labels to ensure exact matches
clean_taxa_names <- function(labels) {
    labels <- gsub("\\.", "_", labels)
    labels <- gsub("-", "_", labels)
    labels <- gsub(" ", "_", labels)
    labels <- gsub("_+", "_", labels) 
    return(labels)
}

calculate_congruence <- function(gene_tree_path, reference_tree) {
    # Load and sanitize gene tree
    gene_tree <- read.tree(gene_tree_path)
    gene_tree$tip.label <- clean_taxa_names(gene_tree$tip.label)
    
    # Ensure unrooted and bifurcating for valid topological comparison
    gene_tree <- multi2di(unroot(gene_tree))
    reference_tree <- multi2di(unroot(reference_tree))
    
    # Identify taxa present in both trees
    common <- intersect(gene_tree$tip.label, reference_tree$tip.label)
    n_taxa <- length(common)
    
    # Minimum 4 taxa needed for unrooted topology comparison
    if(n_taxa < 4) {
        return(c(SD = NA, Score = NA, Taxa = n_taxa))
    }
    
    # Prune trees to the shared taxa set
    p_gene <- keep.tip(gene_tree, common)
    p_ref <- keep.tip(reference_tree, common)
    
    # Calculate Symmetric Difference using ape::dist.topo
    # This is equivalent to Robinson-Foulds distance
    sd_val <- dist.topo(p_gene, p_ref, method = "PH85")
    
    # Normalized Congruence Score (1 - Normalized SD)
    # Max SD distance for unrooted trees is 2 * (n - 3)
    max_sd <- 2 * (n_taxa - 3)
    score <- 1 - (sd_val / max_sd)
    
    return(c(SD = sd_val, Score = score, Taxa = n_taxa))
}

# --- 3. Main Execution ---

cat("🚀 Starting Analysis (ape-only mode)...\n")

if (!file.exists(ref_tree_path)) {
    stop(paste("❌ Reference tree not found at:", ref_tree_path))
}
ref_tree <- read.tree(ref_tree_path)
ref_tree$tip.label <- clean_taxa_names(ref_tree$tip.label)

tree_files <- list.files(path = gene_tree_folder, pattern = "\\.nwk$", full.names = TRUE)
if (length(tree_files) == 0) {
    stop(paste("❌ No .nwk files found in folder:", gene_tree_folder))
}

results <- data.frame(
    Gene = basename(tree_files),
    SD_Distance = NA,
    Congruence_Score = NA,
    Common_Taxa = NA
)

for (i in 1:length(tree_files)) {
    cat(sprintf("[%d/%d] Processing: %s\n", i, length(tree_files), basename(tree_files[i])))
    
    tryCatch({
        stats <- calculate_congruence(tree_files[i], ref_tree)
        results$SD_Distance[i] <- stats["SD"]
        results$Congruence_Score[i] <- stats["Score"]
        results$Common_Taxa[i] <- stats["Taxa"]
    }, error = function(e) {
        cat(sprintf("    ⚠️ Failed to process %s: %s\n", basename(tree_files[i]), conditionMessage(e)))
    })
}

# Sort and Save
results <- results[order(results$Congruence_Score, decreasing = TRUE), ]
write.csv(results, output_csv, row.names = FALSE)

cat("\n✅ Analysis Complete!\n")
print(head(results, 5))
