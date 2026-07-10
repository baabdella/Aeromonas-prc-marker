#!/bin/bash

# Create output folders
mkdir -p ./gene_trees_nwk
mkdir -p ./gene_trees_pdf

echo "Starting Tree Construction and PDF Rendering..."

for f in ./aligned_genes/*_aligned.fasta; do
    base=$(basename "$f" _aligned.fasta)
    
    echo "------------------------------------------"
    echo "Processing Gene: $base"
    
    # 1. Run FastTree (Removed -quote)
    # Use sed to remove any single quotes (') that FastTree might still insert
    fasttree -nt -gtr "$f" | sed "s/'//g" > "./gene_trees_nwk/${base}.nwk"
    
    # 2. Use R to convert Newick to PDF
    Rscript -e "
        library(ape)
        # Check if file is empty before reading
        if (file.size('./gene_trees_nwk/${base}.nwk') > 0) {
            tree <- read.tree('./gene_trees_nwk/${base}.nwk')
            
            # Sanitization inside R to be 100% sure names match
            tree$tip.label <- gsub('[. -]', '_', tree$tip.label)
            
            pdf('./gene_trees_pdf/${base}.pdf', width=10, height=12)
            plot(tree, 
                 main=paste('Phylogeny of', '$base'), 
                 cex=0.8, 
                 edge.width=1.5,
                 no.margin=TRUE)
            add.scale.bar()
            dev.off()
        }
    "
done

echo "Done!"
