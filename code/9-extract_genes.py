from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# --- CONFIGURATION ---
# Path to the folder containing your .fna or .fasta files
genome_folder = "/media/bahaa/Data/Work/2026/validation/raw_data/chromosomes/"
# Path to your coordinate file (Format: ID Start-End)
coord_file = "formatted_coords.txt"
# Name of the output file for your alignment
output_fasta = "for_alignment.fasta"

def extract_sequences():
    # 1. Index all genomes into memory for fast access
    print("Step 1: Indexing genomes from folder...")
    genome_dict = {}
    for filename in os.listdir(genome_folder):
        if filename.endswith((".fna", ".fasta", ".fa")):
            path = os.path.join(genome_folder, filename)
            for record in SeqIO.parse(path, "fasta"):
                # Store by ID (e.g., A_hydrophila_ATCC_7966)
                genome_dict[record.id] = record

    print(f"Loaded {len(genome_dict)} unique sequence IDs.")

    # 2. Process the coordinate list and extract
    print("Step 2: Extracting regions based on coordinates...")
    extracted_records = []
    missing_ids = set()

    with open(coord_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            
            # Expected line: ID Start-End
            parts = line.split()
            if len(parts) < 2:
                print(f"Skipping malformed line: {line}")
                continue
            
            seq_id = parts[0]
            try:
                # Split the Start-End range
                start_str, end_str = parts[1].split("-")
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                print(f"Skipping line with invalid coordinates: {line}")
                continue

            if seq_id in genome_dict:
                parent_seq = genome_dict[seq_id].seq
                
                # Handle Strand and 0-based indexing
                if start < end:
                    # Forward Strand (1-based [start, end] -> 0-based [start-1:end])
                    sub_seq = parent_seq[start-1:end]
                else:
                    # Reverse Strand (Start > End)
                    # We slice the region, then get the reverse complement
                    sub_seq = parent_seq[end-1:start].reverse_complement()
                
                # Create a clean SeqRecord (No description, just the ID)
                new_record = SeqRecord(
                    sub_seq, 
                    id=seq_id, 
                    description=""  # Clearing this prevents the "extracted region..." text
                )
                extracted_records.append(new_record)
            else:
                missing_ids.add(seq_id)

    # 3. Save the results
    if extracted_records:
        SeqIO.write(extracted_records, output_fasta, "fasta")
        print(f"--- SUCCESS ---")
        print(f"Extracted {len(extracted_records)} sequences to: {output_fasta}")
    else:
        print("--- FAILED --- No sequences were extracted.")

    if missing_ids:
        print(f"WARNING: The following {len(missing_ids)} IDs were not found in your genome files:")
        for m_id in sorted(missing_ids):
            print(f"  - {m_id}")

if __name__ == "__main__":
    extract_sequences()
