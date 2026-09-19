import sys
from Bio import SeqIO

INPUT_FASTA = "data/raw_search_hits.fasta" 
OUTPUT_FASTA = "data/radical_SAM_nr.fasta"
EXPECTED_LENGTH = 350                
LENGTH_THRESHOLD = EXPECTED_LENGTH * 0.80  # 80% sequence completeness threshold

def filter_quality_sequences():
    print(f"[INFO] Initiating quality filter (Min length threshold: {LENGTH_THRESHOLD} aa)")
    
    filtered_records = []
    total_count = 0
    passed_count = 0

    try:
        for record in SeqIO.parse(INPUT_FASTA, "fasta"):
            total_count += 1
            
            # keep sequences at least 80% of the expected length
            if len(record.seq) >= LENGTH_THRESHOLD:
                filtered_records.append(record)
                passed_count += 1
                
    except FileNotFoundError:
        print(f"[FATAL] Input file '{INPUT_FASTA}' not found. Terminating.")
        sys.exit(1)

    SeqIO.write(filtered_records, OUTPUT_FASTA, "fasta")
    
    print(f"[INFO] Filtering complete. Evaluated: {total_count} | Retained: {passed_count} | Discarded: {total_count - passed_count}")
    print(f"[INFO] Quality-filtered sequences saved to {OUTPUT_FASTA}")

if __name__ == "__main__":
    filter_quality_sequences()