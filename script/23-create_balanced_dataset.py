import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------
# 1. Paths and constants
# ---------------------------------------------------------------------

master_file = "../data/processed/X_master_dense_prob.csv"
label_file = "../data/processed/y_master.csv"

# Output files
X_out = "../data/processed/X_balanced_dense.csv"
y_out = "../data/processed/y_balanced.csv"

# ---------------------------------------------------------------------
# 2. Load data and drop rows with any NaNs in X
# ---------------------------------------------------------------------
X = pd.read_csv(master_file)
y = pd.read_csv(label_file).squeeze("columns")

assert len(X) == len(y), "Feature matrix and label vector have different lengths!"

nan_mask = X.isna().any(axis=1)
kept_mask = ~nan_mask

X = X.loc[kept_mask].reset_index(drop=True)
y = y.loc[kept_mask].reset_index(drop=True)

print(f"After dropping NaNs: {X.shape[0]:,} rows remain")

# ---------------------------------------------------------------------
# 3. Combine features and true labels
# ---------------------------------------------------------------------
df = X.copy()
df["label"] = y

# ---------------------------------------------------------------------
# 4. Separate positives and negatives
# ---------------------------------------------------------------------
positives = df[df["label"] == 1].copy()
negatives = df[df["label"] == 0].copy()

print(f"Total positives: {len(positives):,}")
print(f"Total true negatives: {len(negatives):,}")

# ---------------------------------------------------------------------
# 5. Sample negatives to balance the dataset
# ---------------------------------------------------------------------
negatives_bal = negatives.sample(
    n=len(positives),
    random_state=42,
    replace=False
).copy()

print(f"Sampled negatives: {len(negatives_bal):,}")

# ---------------------------------------------------------------------
# 6. Combine and shuffle
# ---------------------------------------------------------------------
balanced_df = (
    pd.concat([positives, negatives_bal])
      .sample(frac=1, random_state=42)
      .reset_index(drop=True)
)

X_balanced = balanced_df.drop(columns=["label"])
y_balanced = balanced_df["label"]

# ---------------------------------------------------------------------
# 7. Save
# ---------------------------------------------------------------------
X_balanced.to_csv(X_out, index=False)
y_balanced.to_csv(y_out, index=False)

print("\nSaved:")
print(f"  {X_out}")
print(f"  {y_out}")
print("\n=== Final dataset shape & label counts ===")
print(f"Shape: {X_balanced.shape}")
print(f"Positives: {y_balanced.sum():,}")
print(f"Negatives: {(y_balanced == 0).sum():,}")
