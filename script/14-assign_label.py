# import pandas as pd

# positive_residues = pd.read_csv("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/BioLiP_positives_residue_level_with_duplicates.csv")

# chunk_size = 10000  # Adjust based on your RAM
# chunks = []

# for chunk in pd.read_csv('/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/BioPDB_features.csv', chunksize=chunk_size):
#     processed_chunk = chunk
#     chunks.append(processed_chunk)

# # Combine all chunks if needed
# bioPDF_df = pd.concat(chunks, ignore_index=True)

# # Keep only identifier columns from the positives and drop any duplicates
# pos_keys = (
#     positive_residues[["pdb_id", "chain_id", "residue_number"]]
#     .drop_duplicates()
# )

# # Left-join those keys onto bioPDF_df, tagging matches with label = 1
# bioPDF_df = (
#     bioPDF_df.merge(
#         pos_keys.assign(label=1),
#         on=["pdb_id", "chain_id", "residue_number"],
#         how="left"
#     )
# )


# # 3) Everything that didn’t match gets NaN → set to 0; cast to int
# bioPDF_df["label"] = bioPDF_df["label"].fillna(0).astype(int)

# bioPDF_df.to_csv('/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/BioPDB_features_with_labels.csv', index=False)



import pandas as pd
from pathlib import Path

BIOPDB = Path("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/BioPDB_features.csv")
OUT    = Path("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/BioPDB_features_with_labels.csv")
BIOLIP_POS = Path("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/BioLiP_positives_residue_level_with_duplicates.csv")

# ──────────────────────────────────────────────────────────
# 1) Read only the identifier columns from the positives
#    and build an in-memory lookup that is tiny.
# ──────────────────────────────────────────────────────────
pos_keys = (
    pd.read_csv(BIOLIP_POS,
                usecols=["pdb_id", "chain_id", "residue_number", "renum_residue_number"])
      .drop_duplicates()
)

# Build a Python set for O(1) membership testing
pos_set = set(zip(pos_keys.pdb_id,
                  pos_keys.chain_id,
                  pos_keys.residue_number))

# ──────────────────────────────────────────────────────────
# 2) Stream through BioPDB_features.csv in manageable slices
# ──────────────────────────────────────────────────────────
chunk_size = 50_000          # tune to your RAM
first = True

for chunk in pd.read_csv(BIOPDB, chunksize=chunk_size):
    # vectorised membership test
    keys = list(zip(chunk.pdb_id, chunk.chain_id, chunk.residue_number))
    chunk["label"] = [1 if k in pos_set else 0 for k in keys]

    # write / append to the output file
    chunk.to_csv(OUT, mode="w" if first else "a",
                 index=False, header=first)
    first = False
