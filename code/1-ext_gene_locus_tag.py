import pandas as pd
import sys
import os

def extract_tags_to_folder(input_list_file, csv_file='gene_presence_absence.csv', folder_name='extracted_locus_tags'):
    # 1. Validation
    if not os.path.exists(csv_file):
        print(f"Error: {csv_file} not found.")
        return
    if not os.path.exists(input_list_file):
        print(f"Error: Gene list file '{input_list_file}' not found.")
        return

    # 2. Create the output folder if it doesn't exist
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"Created directory: {folder_name}")

    # 3. Read the gene list
    with open(input_list_file, 'r') as f:
        gene_list = [line.strip() for line in f if line.strip()]

    # 4. Load the Roary CSV
    print(f"Loading {csv_file}...")
    df = pd.read_csv(csv_file, low_memory=False)

    # Define metadata columns to exclude (Strict list)
    metadata_cols = [
        'Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 
        'No. sequences', 'Avg sequences per isolate', 'Genome Fragment', 
        'Order within Fragment', 'Accession', 'Specific Gene Groups', 
        'Average Length (bp)', 'Compound Gene Name'
    ]
    
    isolate_cols = [c for c in df.columns if c not in metadata_cols]

    # 5. Process Genes
    print(f"Starting extraction for {len(gene_list)} genes...")
    
    for target_gene in gene_list:
        gene_row = df[df['Gene'] == target_gene]

        if gene_row.empty:
            print(f"[-] Skipping: '{target_gene}' not found in CSV.")
            continue

        locus_tags = []
        for col in isolate_cols:
            cell_value = gene_row[col].values[0]
            if pd.notna(cell_value):
                # Split handles multiple tags and removes extra whitespace
                tags = str(cell_value).split()
                for t in tags:
                    # Filter out purely numeric entries (metadata error prevention)
                    if not t.replace('.','').isdigit():
                        locus_tags.append(t)

        # 6. Write to file inside the folder
        if locus_tags:
            output_path = os.path.join(folder_name, f"{target_gene}_locus_tags.txt")
            with open(output_path, 'w') as f:
                for tag in locus_tags:
                    f.write(f"{tag}\n")
            print(f"[+] Saved {len(locus_tags)} tags to {output_path}")
        else:
            print(f"[-] No valid locus tags for '{target_gene}'.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <genes_list_file.txt>")
    else:
        extract_tags_to_folder(sys.argv[1])
