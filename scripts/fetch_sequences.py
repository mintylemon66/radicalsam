import sys
from Bio import SeqIO

# --- Configuration ---
IDS_FILE = "radical_SAM_ids.txt"           
DATABASE_FILE = "data/uniprot_trembl.fasta"
OUTPUT_FASTA = "data/raw_search_hits.fasta"    

def fetch_local_sequences():
    print(f"[INFO] Initializing sequence extraction from {DATABASE_FILE}")
    
    try:
        with open(IDS_FILE, "r") as f:
            # Hash map implementation for O(1) lookup times
            target_ids = set(line.strip() for line in f if line.strip())
    except FileNotFoundError:
        print(f"[FATAL] ID list '{IDS_FILE}' not found. Terminating.")
        sys.exit(1)

    print(f"[INFO] Loaded {len(target_ids)} target IDs.")
    
    found_records = []
    
    try:
        # Iterate through local database and extract matching records
        for record in SeqIO.parse(DATABASE_FILE, "fasta"):
            if record.id in target_ids:
                found_records.append(record)
    except FileNotFoundError:
         print(f"[FATAL] Database file '{DATABASE_FILE}' not found. Terminating.")
         sys.exit(1)
            
    print(f"[INFO] Extraction complete. {len(found_records)} matching sequences retrieved.")
    
    SeqIO.write(found_records, OUTPUT_FASTA, "fasta")
    print(f"[INFO] Sequences written to {OUTPUT_FASTA}")

if __name__ == "__main__":
    fetch_local_sequences()