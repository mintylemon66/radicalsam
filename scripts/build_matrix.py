import sys
import pandas as pd

# --- Configuration ---
METADATA_FILE = "data/metadata_org.tsv"
NEIGHBOR_FILE = "data/neighborhood_genes.tsv"
ANNOTATION_FILE = "data/neighborhoods_annotated.tsv"
OUTPUT_MATRIX = "data/neighborhood_matrix.tsv"

def main():
    print("--- Roadmap B, Step 4: Building the Master Matrix ---")

    # 1. Load the raw neighborhood mappings
    try:
        print("[INFO] Loading neighborhood maps...")
        df_neighbors = pd.read_csv(NEIGHBOR_FILE, sep="\t")
    except FileNotFoundError:
        print(f"[FATAL] Missing {NEIGHBOR_FILE}. Please wait for your overnight script to finish.")
        sys.exit(1)

    # 2. Unpivot the comma-separated genes into individual structural rows
    print("[INFO] Exploding gene lists into structural matrix...")
    matrix_rows = []
    
    for index, row in df_neighbors.iterrows():
        rs_id = row['Radical_SAM_ID']
        contig = row['Genomic_Contig_ID']
        genes_str = str(row['Neighbor_Genes'])
        
        if genes_str and genes_str not in ["nan", "None", "[No_Genes_Found]", "none_found"]:
            gene_list = [g.strip() for g in genes_str.split(",") if g.strip()]
            
            # Find the center index to calculate relative genomic distance rank (-5 to +5)
            center_idx = len(gene_list) // 2 
            
            for i, gene in enumerate(gene_list):
                rank = i - center_idx
                matrix_rows.append({
                    "Sequence_ID": rs_id,  # Named Sequence_ID to match your metadata file
                    "Genomic_Contig_ID": contig,
                    "Neighbor_Rank": rank,
                    "Gene_ID": gene
                })

    df_matrix = pd.DataFrame(matrix_rows)
    print(f"[INFO] Generated base matrix with {len(df_matrix)} gene rows.")

    # 3. Merge with Annotations (InterProScan data)
    try:
        print("[INFO] Merging functional annotations...")
        # InterProScan TSV files have no headers by default. We assign names to columns we need.
        # Column 0: Protein/Gene ID, Column 5: Signature Description, Column 11/12: InterPro ID and Desc
        df_anno = pd.read_csv(ANNOTATION_FILE, sep="\t", header=None, usecols=[0, 5, 11, 12], 
                            names=['Gene_ID', 'Signature_Desc', 'InterPro_ID', 'InterPro_Desc'])
        
        # Deduplicate in case InterPro found multiple domains within a single neighbor gene
        df_anno = df_anno.drop_duplicates(subset=['Gene_ID'])
        
        df_matrix = df_matrix.merge(df_anno, on="Gene_ID", how="left")
    except FileNotFoundError:
        print(f"[WARNING] {ANNOTATION_FILE} not found. Skipping functional merge for now.")

    # 4. Merge with Organism Metadata using exact user headers
    try:
        print("[INFO] Merging organism metadata...")
        df_meta = pd.read_csv(METADATA_FILE, sep="\t")
        
        # Merge directly on your "Sequence_ID" column
        df_matrix = df_matrix.merge(df_meta, on="Sequence_ID", how="left")
    except FileNotFoundError:
        print(f"[WARNING] {METADATA_FILE} not found. Skipping metadata merge.")

    # 5. Save the complete master dataset
    df_matrix.to_csv(OUTPUT_MATRIX, sep="\t", index=False)
    print(f"\n[SUCCESS] Master matrix successfully saved to {OUTPUT_MATRIX}!")
    print(f"[SUCCESS] Final Dimensions: {df_matrix.shape[0]} rows x {df_matrix.shape[1]} columns.")

if __name__ == "__main__":
    main()