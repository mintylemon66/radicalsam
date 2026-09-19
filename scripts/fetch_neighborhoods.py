import sys
import time
import pandas as pd
from Bio import Entrez
import urllib.error

Entrez.email = "mintyl@bu.edu"  # NCBI requires an email to use their API
INPUT_FILE = "data/raw_neighborhoods.txt"
OUTPUT_FILE = "data/neighborhood_genes.tsv"
BATCH_SIZE = 100  # contigs per Entrez request

import http.client

def fetch_ncbi_features_batch(contig_ids, max_retries=3):
    """Fetch feature tables, retrying when the connection drops."""
    id_string = ",".join(contig_ids)
    
    for attempt in range(max_retries):
        try:
            handle = Entrez.efetch(db="nucleotide", id=id_string, rettype="ft", retmode="text")
            raw_data = handle.read()
            handle.close()
            return raw_data
            
        except http.client.IncompleteRead as e:
            print(f"  [WARNING] NCBI hung up (File too large). Retrying {attempt + 1}/{max_retries}...")
            time.sleep(3)
        except urllib.error.HTTPError as e:
            print(f"  [ERROR] NCBI API rejected batch: {e.code}. Retrying...")
            time.sleep(3)
        except Exception as e:
            print(f"  [ERROR] Batch connection failed: {e}. Retrying...")
            time.sleep(3)
            
    print("  [FATAL] Batch failed completely after 3 retries. Skipping chunk.")
    return ""

def parse_feature_tables(raw_text):
    """Parse NCBI multi-record feature table text."""
    if not raw_text.strip():
        return {}
        
    records = raw_text.split(">Feature ")
    contig_to_genes = {}
    
    for record in records:
        if not record.strip():
            continue
        lines = record.splitlines()
        
        # strip gb| / ref| prefixes and trailing pipes from the header
        raw_id = lines[0].split()[0]
        if "|" in raw_id:
            parts = raw_id.split("|")
            # header is usually db|accession|
            raw_id = parts[1] if len(parts) > 1 else parts[0]
            
        clean_first_line = raw_id.split(".")[0]  # drop the version suffix
        
        # protein IDs, tolerant of tab vs space differences
        genes = []
        for line in lines:
            if "protein_id" in line:
                parts = line.strip().split()
                if len(parts) >= 2 and parts[0] == "protein_id":
                    genes.append(parts[1])
                    
        contig_to_genes[clean_first_line] = genes
        
    return contig_to_genes

def main():
    print("--- Roadmap B, Step 1: Bulk Batch Neighborhood Extraction ---")
    
    try:
        df = pd.read_csv(INPUT_FILE, sep="\t")
        print(f"[INFO] Loaded {len(df)} genomic contigs.")
    except FileNotFoundError:
        print(f"[FATAL] Could not find {INPUT_FILE}.")
        sys.exit(1)

    valid_contigs = df.dropna(subset=['Genomic_Contig_ID']).copy()
    total_records = len(valid_contigs)
    print(f"[INFO] Found {total_records} valid contigs processing...")
    
    # strip version numbers so the lookups match
    valid_contigs['Clean_Contig'] = valid_contigs['Genomic_Contig_ID'].str.strip().str.split('.').str[0]
    
    with open(OUTPUT_FILE, "w") as out_file:
        out_file.write("Radical_SAM_ID\tGenomic_Contig_ID\tNeighbor_Genes\n")
        
        for start_idx in range(0, total_records, BATCH_SIZE):
            chunk = valid_contigs.iloc[start_idx : start_idx + BATCH_SIZE]
            contig_list = chunk['Clean_Contig'].tolist()
            
            print(f"[INFO] Processing records {start_idx} to {start_idx + len(chunk)} of {total_records}...")
            
            raw_features = fetch_ncbi_features_batch(contig_list)
            
            # if a whole batch fails, keep what parsed and mark the rest empty
            parsed_neighborhoods = parse_feature_tables(raw_features) if raw_features else {}
            
            # write each batch as soon as it finishes
            for index, row in chunk.iterrows():
                clean_c = row['Clean_Contig']
                orig_c = row['Genomic_Contig_ID']
                rs_id = row['Clean_Accession']
                
                genes = parsed_neighborhoods.get(clean_c, [])
                gene_str = ",".join(genes) if genes else "[No_Genes_Found]"
                
                out_file.write(f"{rs_id}\t{orig_c}\t{gene_str}\n")
            
            # NCBI allows 3 requests per second
            time.sleep(0.4)
            
    print(f"\n[SUCCESS] Bulk extraction complete! All entries processed in {OUTPUT_FILE}")

if __name__ == "__main__":
    main()