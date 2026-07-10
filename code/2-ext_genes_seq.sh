# 1. Create the output folder first
mkdir -p ./extracted_genes

# 2. Loop through your gene list
for x in $(cat my_genes.txt)
do
    echo "Processing $x..."
    
    # Check if the locus tag file exists before running seqkit
    if [ -f "./extracted_locus_tags/${x}_locus_tags.txt" ]; then
        
        seqkit grep -f "./extracted_locus_tags/${x}_locus_tags.txt" \
            ../1_bakta_result/*/*.ffn \
            -o "./extracted_genes/${x}_genes.fa"
            
    else
        echo "Warning: ./extracted_locus_tags/${x}_locus_tags.txt not found."
    fi
done
