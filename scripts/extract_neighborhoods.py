import sys
import time
import json
import urllib.request
import urllib.parse
from Bio import SeqIO

# --- Configuration ---
INPUT_FASTA = "data/radical_SAM_nr_derep.fasta"
OUTPUT_TXT = "data/raw_neighborhoods.txt"
BATCH_SIZE = 50  # 50 is the mathematical sweet spot to avoid URL character limits

def clean_protein_id(raw_id):
    """Extracts the core accession from a UniProt/TrEMBL formatted header."""
    if "|" in raw_id:
        return raw_id.split("|")[1]
    return raw_id

def fetch_batch_from_uniprot(batch_ids):
    """Queries the UniProt API for up to 50 proteins simultaneously."""
    # Construct a search query: accession:ID1 OR accession:ID2 ...
    query_str = " OR ".join([f"accession:{acc}" for acc in batch_ids])
    encoded_query = urllib.parse.quote(query_str)
    
    # We use the search endpoint to bypass the 1-by-1 lookup
    url = f"https://rest.uniprot.org/uniprotkb/search?query={encoded_query}&format=json&size={BATCH_SIZE}"
    
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode('utf-8'))
            
        # Parse the JSON payload for all proteins in the batch at once
        results = {}
        for entry in data.get("results", []):
            acc = entry["primaryAccession"]
            contig = None
            cross_references = entry.get("uniProtKBCrossReferences", [])
            for xref in cross_references:
                if xref["database"] == "EMBL":
                    contig = xref["id"]
                    break
            results[acc] = contig
            
        return results

    except Exception as e:
        print(f"  [ERROR] Batch connection failed: {e}")
        return {}

def extract_neighborhoods():
    print(f"[INFO] Initiating High-Speed Batch Genomic Extraction")
    
    try:
        records = list(SeqIO.parse(INPUT_FASTA, "fasta"))
        print(f"[INFO] Loaded {len(records)} Radical SAM candidates.")
    except FileNotFoundError:
        print(f"[FATAL] Input file '{INPUT_FASTA}' not found. Terminating.")
        sys.exit(1)

    # --- THE AUTO-RESUME FEATURE ---
    # Automatically counts how many you have already done so it doesn't start over!
    try:
        with open(OUTPUT_TXT, "r") as f:
            # Subtract 1 to account for the header row
            completed_lines = max(0, len(f.readlines()) - 1)
    except FileNotFoundError:
        completed_lines = 0

    print(f"[INFO] Detected {completed_lines} already processed sequences.")
    print(f"[INFO] Resuming extraction from sequence {completed_lines}...")
    
    records_to_process = records[completed_lines:]
    successful_mappings = 0

    # Open in APPEND mode to safely protect your existing data
    with open(OUTPUT_TXT, "a") as out_file:
        # If it is a completely fresh run, write the header
        if completed_lines == 0:
            out_file.write("Radical_SAM_ID\tClean_Accession\tGenomic_Contig_ID\tNeighbor_Genes\n")
        
        # Loop through the remaining sequences in chunks of 50
        for i in range(0, len(records_to_process), BATCH_SIZE):
            batch_records = records_to_process[i:i + BATCH_SIZE]
            
            # Map the clean ID back to the messy FASTA ID
            id_map = {clean_protein_id(rec.id): rec.id for rec in batch_records}
            clean_ids = list(id_map.keys())
            
            print(f"[INFO] Processing Batch: Sequences {completed_lines + i} to {completed_lines + i + len(batch_records)}...")
            batch_results = fetch_batch_from_uniprot(clean_ids)
            
            # Write successful hits to the file
            for clean_id, contig_id in batch_results.items():
                if contig_id:
                    original_id = id_map[clean_id]
                    out_file.write(f"{original_id}\t{clean_id}\t{contig_id}\t[Pending_Neighborhood_Extraction]\n")
                    successful_mappings += 1
            
            # Pause for 1 second to respect UniProt API limits
            time.sleep(1)

    print(f"\n[SUCCESS] Batch extraction complete!")

if __name__ == "__main__":
    extract_neighborhoods()