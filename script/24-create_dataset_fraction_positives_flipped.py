import pandas as pd
from pathlib import Path
import numpy as np

# ---------------------------------------------------------------------
# 1. Paths and constants
# ---------------------------------------------------------------------
DATA_DIR = Path("/home/mpradhan/Intern_Research_Project/data/")
master_file = DATA_DIR / "X_master_dense_prob.csv"
label_file = DATA_DIR / "y_master.csv"

# Output files
X_out      = DATA_DIR / "X_balanced_validation.csv"
y_true_out = DATA_DIR / "y_balanced_true_validation.csv"
y_pu_out   = DATA_DIR / "y_balanced_validation.csv"

# How many positives to flip to unlabeled
FLIP_FRACTION = 0.10    # 10 % of positives will become unlabeled
N_UNLABELED   = 836_000 # How many unlabeled (true negatives) to keep

# ---------------------------------------------------------------------
# 2. Load data and drop rows with any NaNs in X
# ---------------------------------------------------------------------
X = pd.read_csv(master_file)
y = pd.read_csv(label_file).squeeze("columns")

assert len(X) == len(y), "Feature matrix and label vector have different lengths!"

nan_mask   = X.isna().any(axis=1)
kept_mask  = ~nan_mask

X = X.loc[kept_mask].reset_index(drop=True)
y = y.loc[kept_mask].reset_index(drop=True)

print(f"After dropping NaNs: {X.shape[0]:,} rows remain")

# ---------------------------------------------------------------------
# 3. Combine features and true labels
# ---------------------------------------------------------------------
df = X.copy()
df["true_label"] = y

# ---------------------------------------------------------------------
# 4. Separate true positives and true negatives
# ---------------------------------------------------------------------
positives = df[df["true_label"] == 1].copy()
negatives = df[df["true_label"] == 0].copy()

print(f"Total positives: {len(positives):,}")
print(f"Total true negatives: {len(negatives):,}")

# ---------------------------------------------------------------------
# 5. Flip 10 % of positives to unlabeled for PU learning
# ---------------------------------------------------------------------
n_flip       = int(len(positives) * FLIP_FRACTION)
flipped_idx  = positives.sample(n=n_flip, random_state=42).index

positives["pu_label"] = 1
positives.loc[flipped_idx, "pu_label"] = 0

# ---------------------------------------------------------------------
# 6. RANDOMLY sample true negatives to serve as unlabeled
# ---------------------------------------------------------------------
negatives_sub = negatives.sample(
    n=N_UNLABELED,
    random_state=42,   # for reproducibility; adjust or remove for a new draw
    replace=False
).copy()
negatives_sub["pu_label"] = 0  # always unlabeled

print(f"Randomly selected unlabeled: {len(negatives_sub):,}")

# ---------------------------------------------------------------------
# 7. Combine positives + sampled negatives, then shuffle
# ---------------------------------------------------------------------
new_df = (
    pd.concat([positives, negatives_sub])
      .sample(frac=1, random_state=42)
      .reset_index(drop=True)
)

# Split out features and labels
X_new  = new_df.drop(columns=["true_label", "pu_label"])
y_true = new_df["true_label"]   # 1 = true positive, 0 = true negative
y_pu   = new_df["pu_label"]     # 1 = observed positive, 0 = unlabeled

# ---------------------------------------------------------------------
# 8. Save results
# ---------------------------------------------------------------------
X_new.to_csv(X_out,      index=False)
y_true.to_csv(y_true_out, index=False)
y_pu.to_csv(y_pu_out,     index=False)

print("\nSaved:")
print(f"  {X_out}")
print(f"  {y_true_out}")
print(f"  {y_pu_out}")

print("\n=== Final dataset shape & label counts ===")
print(f"Shape: {X_new.shape}")
print(f"True Positives: {y_true.sum():,}")
print(f"PU Positives:  {y_pu.sum():,}")
print(f"PU Unlabeled: {(y_pu == 0).sum():,}")

# Detailed breakdown
true_pos_pu_pos  = ((y_true == 1) & (y_pu == 1)).sum()
true_pos_pu_unlb = ((y_true == 1) & (y_pu == 0)).sum()
true_neg_pu_unlb = ((y_true == 0) & (y_pu == 0)).sum()

print("\n=== Detailed breakdown ===")
print(f" True positives (y_true == 1): {y_true.sum():,}")
print(f"   ├─ still labeled positive: {true_pos_pu_pos:,}")
print(f"   └─ flipped to unlabeled:   {true_pos_pu_unlb:,}")
print(f" True negatives (y_true == 0): {(y_true == 0).sum():,}")
print(f"   └─ remain unlabeled:       {true_neg_pu_unlb:,}")
print(f" Final PU labels → Positives: {y_pu.sum():,}, Unlabeled: {(y_pu == 0).sum():,}")
