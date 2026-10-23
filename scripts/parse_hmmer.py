import os
import sys

# --- Configuration ---
INPUT_TSV = "data/radical_SAM_all_hits.txt"
OUTPUT_IDS_FILE = "data/radical_SAM_ids.txt"
E_VALUE_CUTOFF = 1e-10

def parse_hmmer_tblout():
    print(f"[INFO] Parsing HMMER table output: {INPUT_TSV}")
    
    if not os.path.exists(INPUT_TSV):
        print(f"[FATAL] File '{INPUT_TSV}' not found. Terminating.")
        sys.exit(1)

    # Use a set for automatic deduplication
    target_ids = set()
    
    with open(INPUT_TSV, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
                
            parts = line.split()
            if len(parts) < 13:
                continue
                
            # HMMER tblout format: column 0 = target name, column 4 = full sequence E-value
            target_name = parts[0]
            try:
                e_value = float(parts[4])
            except ValueError:
                continue
                
            # Isolate sequences meeting the strict confidence threshold
            if e_value <= E_VALUE_CUTOFF:
                target_ids.add(target_name)

    with open(OUTPUT_IDS_FILE, "w") as out_f:
        for target_id in target_ids:
            out_f.write(f"{target_id}\n")
            
    print(f"[INFO] Parsing complete. {len(target_ids)} significant sequences identified.")
    print(f"[INFO] Sequence IDs written to {OUTPUT_IDS_FILE}")

if __name__ == "__main__":
    parse_hmmer_tblout()