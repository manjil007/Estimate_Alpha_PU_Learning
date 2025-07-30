# prep_split_pdb_ids.ipynb
import pandas as pd
import numpy as np
from pathlib import Path

# ------------------------------------------------------------------
# 1. Load the two BioLiP residue‑level tables
# ------------------------------------------------------------------
pos_df  = pd.read_csv("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/BioLiP_positives_residue_level_with_duplicates.csv")
unlab_df = pd.read_csv("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/BioLiP_unlabeled_residues_with_duplicates.csv")

# ------------------------------------------------------------------
# 2. Collect **all** unique PDB IDs (positives + unlabeled)
# ------------------------------------------------------------------
all_pdb_ids = pd.concat([pos_df["pdb_id"], unlab_df["pdb_id"]]).drop_duplicates().sort_values().reset_index(drop=True)

print(f"Total unique PDB IDs: {len(all_pdb_ids)}")

# ------------------------------------------------------------------
# 3. Split into three roughly equal parts
# ------------------------------------------------------------------
parts = np.array_split(all_pdb_ids, 3)

out_dir = Path("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/raw")
out_dir.mkdir(exist_ok=True, parents=True)

for i, part in enumerate(parts, 1):
    outfile = out_dir / f"pdb_ids_part{i}.csv"
    part.to_csv(outfile, index=False, header=["pdb_id"])
    print(f"[✓] Saved {len(part)} IDs → {outfile}")
