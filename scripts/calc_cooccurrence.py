import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import sys

# --- Configuration ---
INPUT_FILE = "data/raw_neighborhoods.txt"
OUTPUT_MATRIX = "data/jaccard_matrix.csv"
TOP_N_GENES = 200 # Filter to the most frequent neighbors to save memory

def load_and_clean_data():
    """Loads the neighborhood data and extracts the neighbor gene lists."""
    print(f"[INFO] Loading data from {INPUT_FILE}...")
    try:
        # Read the tab-separated file
        df = pd.read_csv(INPUT_FILE, sep="\t")
        
        # Drop rows where neighborhood extraction is pending or failed
        df = df.dropna(subset=['Neighbor_Genes'])
        df = df[df['Neighbor_Genes'] != "[Pending_Neighborhood_Extraction]"]
        
        print(f"[SUCCESS] Loaded {len(df)} fully extracted neighborhoods.")
        return df
    except FileNotFoundError:
        print(f"[FATAL] Could not find {INPUT_FILE}.")
        sys.exit(1)

def build_binary_matrix(df):
    """Converts lists of neighbor genes into a binary presence/absence matrix."""
    print("[INFO] Building binary occurrence matrix...")
    
    # Explode the comma-separated strings of genes into a list, then into separate rows
    # Assuming 'Neighbor_Genes' looks like "GeneA,GeneB,GeneC"
    df['Neighbor_Genes'] = df['Neighbor_Genes'].str.split(',')
    exploded = df.explode('Neighbor_Genes')
    
    # Get the top N most frequent neighbor genes to prevent memory crash
    top_genes = exploded['Neighbor_Genes'].value_counts().head(TOP_N_GENES).index
    exploded = exploded[exploded['Neighbor_Genes'].isin(top_genes)]
    
    # Create a crosstab (binary matrix: Rows = Radical SAM Contigs, Cols = Neighbor Genes)
    binary_matrix = pd.crosstab(exploded['Genomic_Contig_ID'], exploded['Neighbor_Genes'])
    
    # Convert counts to binary (1 if present, 0 if absent)
    binary_matrix = (binary_matrix > 0).astype(int)
    
    print(f"[INFO] Matrix shape: {binary_matrix.shape} (Contigs x Top {TOP_N_GENES} Genes)")
    return binary_matrix

def calculate_jaccard(binary_matrix):
    """Calculates the pairwise Jaccard similarity between neighbor genes."""
    print("[INFO] Calculating Jaccard similarity...")
    
    # Transpose the matrix so rows are Genes and columns are Contigs
    gene_matrix = binary_matrix.T
    
    # pdist calculates pairwise distances. Jaccard distance is 1 - Jaccard similarity.
    # We subtract from 1 to get the Similarity Index.
    jaccard_distances = pdist(gene_matrix.values, metric='jaccard')
    jaccard_similarity = 1 - squareform(jaccard_distances)
    
    # Convert back to a beautifully labeled Pandas DataFrame
    similarity_df = pd.DataFrame(
        jaccard_similarity, 
        index=gene_matrix.index, 
        columns=gene_matrix.index
    )
    
    return similarity_df

def main():
    print("--- Phase 2: Co-occurrence Matrix Calculation ---")
    
    df = load_and_clean_data()
    binary_matrix = build_binary_matrix(df)
    similarity_df = calculate_jaccard(binary_matrix)
    
    # Save the final matrix
    print(f"[INFO] Saving matrix to {OUTPUT_MATRIX}...")
    similarity_df.to_csv(OUTPUT_MATRIX)
    print("[SUCCESS] Pipeline complete! Your evolutionary partners are ready for review.")

if __name__ == "__main__":
    main()