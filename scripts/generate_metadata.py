import time
import re
import csv
import sys
from Bio import Entrez
from Bio import SeqIO

# --- Configuration ---
Entrez.email = "mintyl@bu.edu" 
INPUT_FASTA = "data/radical_SAM_nr_derep.fasta" 
OUTPUT_TSV = "data/metadata_org.tsv"

def extract_tax_id(fasta_header):
    """Extracts NCBI Taxonomy ID from UniProt OX field."""
    match = re.search(r'OX=(\d+)', fasta_header)
    return match.group(1) if match else None

def get_taxonomy_info(tax_id):
    """Retrieves Scientific Name and Phylum from NCBI Taxonomy API."""
    try:
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        if records:
            record = records[0]
            scientific_name = record.get("ScientificName", "Unknown")
            
            # Parse lineage for Phylum rank
            phylum = "Unknown"
            for clade in record.get("LineageEx", []):
                if clade.get("Rank") == "phylum":
                    phylum = clade.get("ScientificName")
                    break

            return scientific_name, phylum
            
    except Exception as e:
        print(f"[ERROR] API request failed for TaxID {tax_id}: {e}")

    return "Unknown", "Unknown"

def main():
    print(f"[INFO] Initiating metadata extraction for {INPUT_FASTA}")
    
    # Cache to prevent redundant API calls for highly represented clades
    tax_cache = {}

    try:
        with open(OUTPUT_TSV, "w", newline='') as tsv_file:
            writer = csv.writer(tsv_file, delimiter='\t')
            writer.writerow(["Sequence_ID", "Organism_Name", "Phylum", "Habitat", "GC_Percent"])

            for record in SeqIO.parse(INPUT_FASTA, "fasta"):
                tax_id = extract_tax_id(record.description)

                if tax_id:
                    if tax_id not in tax_cache:
                        time.sleep(0.35) # NCBI API rate limit compliance
                        tax_cache[tax_id] = get_taxonomy_info(tax_id)

                    name, phylum = tax_cache[tax_id]
                else:
                    name, phylum = "Not Found", "Not Found"

                # Habitat and GC% reserved for post-processing curation
                writer.writerow([record.id, name, phylum, "TODO", "TODO"])
                
        print(f"[INFO] Metadata table successfully written to {OUTPUT_TSV}")

    except FileNotFoundError:
        print(f"[FATAL] Input file {INPUT_FASTA} not found. Terminating.")
        sys.exit(1)

if __name__ == "__main__":
    main()