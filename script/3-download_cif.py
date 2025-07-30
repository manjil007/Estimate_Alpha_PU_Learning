import os
import requests
import pandas as pd
from pathlib import Path
from tqdm import tqdm

def download_cif_files(residue_dfs, output_dir):
    """
    Download mmCIF files for all unique pdb_ids in the provided residue dataframes.

    Args:
        residue_dfs (list of pd.DataFrame): List of dataframes containing 'pdb_id' columns.
        output_dir (str): Path to store downloaded .cif files.
    """
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Collect unique pdb_ids from all dataframes
    pdb_ids = set()
    for df in residue_dfs:
        pdb_ids.update(df["pdb_id"].dropna().str.lower().unique())

    print(f"Preparing to download {len(pdb_ids)} .cif files...")

    failed = []
    for pdb_id in tqdm(pdb_ids):
        cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        cif_path = os.path.join(output_dir, f"{pdb_id}.cif")

        if os.path.exists(cif_path):
            continue  # Skip if already downloaded

        try:
            response = requests.get(cif_url, timeout=10)
            response.raise_for_status()
            with open(cif_path, 'wb') as f:
                f.write(response.content)
        except Exception as e:
            print(f"Failed to download {pdb_id}: {e}")
            failed.append(pdb_id)

    print(f"\nDownloaded {(len(pdb_ids) - len(failed))} files successfully.")
    if failed:
        print(f"Failed to download {len(failed)} files: {failed}")

# Example usage
if __name__ == "__main__":
    # Replace with actual file paths if needed
    residue_level_df = pd.read_csv("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/BioLiP_positives_residue_level_with_duplicates.csv")
    unlabeled_residue_df = pd.read_csv("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/BioLiP_unlabeled_residues_with_duplicates.csv")

    download_cif_files([residue_level_df, unlabeled_residue_df], output_dir="/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/raw/structures_cif")
