mkdir -p ./aligned_genes

# Runs mafft on 4 cores simultaneously
ls ./renamed_genes/*.fa | parallel -j 4 "mafft --auto {} > ./aligned_genes/{/.}_aligned.fasta"
