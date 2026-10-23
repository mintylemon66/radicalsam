import sys
import time
import urllib.error
from Bio import Entrez

# --- Configuration ---
Entrez.email = "mintyl@bu.edu"  # Use your actual registered email
INPUT_FILE = "data/neighborhood_genes.tsv"
OUTPUT_FASTA = "data/unique_neighbors.fasta"
BATCH_SIZE = 200  # Safe batch size for fetching raw FASTA sequences

def main():
    print("--- Roadmap B, Step 2: Extracting Unique Sequences ---")
    
    # 1. Collect all unique gene IDs from your neighborhood map
    unique_ids = set()
    try:
        print(f"[INFO] Reading {INPUT_FILE}...")
        with open(INPUT_FILE, "r") as f:
            header = f.readline()  # Skip header row
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                
                # Grab the comma-separated list of neighbor genes
                genes_str = parts[2]
                if genes_str and genes_str not in ["nan", "None", "[No_Genes_Found]", "none_found"]:
                    # Split them up, clean white spaces, and add to our set
                    for gene in genes_str.split(","):
                        clean_gene = gene.strip()
                        if clean_gene:
                            unique_ids.add(clean_gene)
    except FileNotFoundError:
        print(f"[FATAL] Could not find {INPUT_FILE}. Wait for your overnight script to finish!")
        sys.exit(1)
        
    id_list = list(unique_ids)
    total_ids = len(id_list)
    print(f"[SUCCESS] Found {total_ids} unique neighbor protein IDs to annotate.")
    
    if total_ids == 0:
        print("[WARNING] No neighbor IDs found. Exiting.")
        sys.exit(0)

    # 2. Batch download their FASTA sequences from NCBI Protein Database
    print(f"[INFO] Downloading FASTA sequences into {OUTPUT_FASTA}...")
    
    with open(OUTPUT_FASTA, "w") as out_file:
        for i in range(0, total_ids, BATCH_SIZE):
            chunk = id_list[i:i + BATCH_SIZE]
            print(f"  -> Fetching sequences {i} to {min(i + BATCH_SIZE, total_ids)}...")
            
            id_string = ",".join(chunk)
            max_retries = 3
            
            for attempt in range(max_retries):
                try:
                    # Query the PROTEIN database for the FASTA formats
                    handle = Entrez.efetch(db="protein", id=id_string, rettype="fasta", retmode="text")
                    fasta_data = handle.read()
                    handle.close()
                    
                    # Write directly to our master FASTA file
                    out_file.write(fasta_data)
                    break  # Success, break the retry loop
                    
                except Exception as e:
                    print(f"    [WARNING] Connection failed ({e}). Retrying loop {attempt + 1}/{max_retries}...")
                    time.sleep(3)
            else:
                print(f"    [ERROR] Failed to download chunk starting at index {i}. Skipping.")
            
            # Polite pause to stay perfectly under the NCBI rate-limiting firewalls
            time.sleep(0.4)

    print(f"\n[SUCCESS] All unique sequences compiled into {OUTPUT_FASTA}!")

if __name__ == "__main__":
    main()