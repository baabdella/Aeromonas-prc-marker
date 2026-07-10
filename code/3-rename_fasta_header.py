import pandas as pd
import os
import sys
import glob

def rename_genes_clean(fasta_file, csv_file='gene_presence_absence.csv', output_folder='renamed_genes'):
    if not os.path.exists(csv_file) or not os.path.exists(fasta_file):
        return

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    base_name = os.path.basename(fasta_file)
    output_file = os.path.join(output_folder, base_name.replace(".fa", "_renamed.fa").replace(".fasta", "_renamed.fasta"))
    
    # 1. Load Roary CSV
    df = pd.read_csv(csv_file, low_memory=False)
    
    # 2. Identify Isolate Columns
    metadata_cols = [
        'Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 
        'No. sequences', 'Avg sequences per isolate', 'Genome Fragment', 
        'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 
        'QC', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc'
    ]
    strain_cols = [c for c in df.columns if c not in metadata_cols]

    # 3. Create Mapping: LocusTag -> StrainName
    tag_to_strain = {}
    for strain in strain_cols:
        col_data = df[strain].dropna().astype(str)
        for entry in col_data:
            for item in entry.replace('\t', ' ').split():
                tag_to_strain[item.strip()] = strain

    # 4. Process the FASTA
    count = 0
    with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                # Get the Locus Tag from the header
                original_id = line.strip().split()[0][1:]
                
                # Direct Lookup
                strain_name = tag_to_strain.get(original_id)
                
                if strain_name:
                    # STRICT SANITIZATION
                    # We replace dots/spaces/hyphens with underscores 
                    # This prevents names from being "cut off" in Newick files
                    clean_name = strain_name.replace(" ", "_").replace(".", "_").replace("-", "_")
                    f_out.write(f">{clean_name}\n")
                    count += 1
                else:
                    # Keep original if not found in Roary
                    f_out.write(f">{original_id}\n")
            else:
                f_out.write(line)

    print(f"File: {base_name} | Renamed {count} sequences.")

if __name__ == "__main__":
    input_path = sys.argv[1] if len(sys.argv) > 1 else "extracted_genes"
    if os.path.isdir(input_path):
        fasta_files = glob.glob(os.path.join(input_path, "*.fa*"))
        for f in fasta_files:
            rename_genes_clean(f)
    else:
        rename_genes_clean(input_path)
