from Bio import SeqIO
import math

INPUT_FILE = "data/clean_neighbors.fasta"

def main():
    print("[INFO] Reading the massive 223MB FASTA file...")
    records = list(SeqIO.parse(INPUT_FILE, "fasta"))
    
    total_seqs = len(records)
    chunk_size = math.ceil(total_seqs / 4) # Split into 4 even pieces
    
    print(f"[INFO] Found {total_seqs} total sequences. Splitting into chunks of {chunk_size}...")
    
    for i in range(4):
        start_idx = i * chunk_size
        end_idx = start_idx + chunk_size
        chunk = records[start_idx:end_idx]
        
        output_name = f"data/clean_neighbors_part_{i+1}.fasta"
        SeqIO.write(chunk, output_name, "fasta")
        print(f"[SUCCESS] Saved {len(chunk)} sequences to {output_name}")

if __name__ == "__main__":
    main()
    