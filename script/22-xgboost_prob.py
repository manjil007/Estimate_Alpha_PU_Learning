#!/usr/bin/env python
# make_probs.py
"""
Train an XGBoost model on the master data (positive vs. unlabeled)
and write a probability column back to disk as X_master_dense_prob.csv.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from xgboost import XGBClassifier

# ----------------------------------------------------------------------
# 1.  Load data
# ----------------------------------------------------------------------


X_path = "../data/processed/X_master_dense.csv"
y_path = "../data/processed/y_master.csv"
out_path = "../data/processed/X_master_dense_prob.csv"

print(f"Loading\n  X → {X_path}\n  y → {y_path}")

# X: rows must be in exactly the same order as y
X = pd.read_csv(X_path)
y = pd.read_csv(y_path).squeeze("columns")              # Series

assert len(X) == len(y), "X and y row counts differ!"
n_pos = int((y == 1).sum())
n_neg = int((y == 0).sum())

print(f"Data shape: {X.shape}  (positives={n_pos}, unlabeled={n_neg})")

# ----------------------------------------------------------------------
# 2.  Fit XGBoost
# ----------------------------------------------------------------------
scale_pos_weight = n_neg / max(n_pos, 1)                # avoid div-by-zero
model = XGBClassifier(
    objective="binary:logistic",
    eval_metric="logloss",
    n_estimators=600,
    learning_rate=0.05,
    max_depth=4,
    subsample=0.8,
    colsample_bytree=0.8,
    scale_pos_weight=scale_pos_weight,
    n_jobs=-1,
    random_state=42,
)

print("Training model …")
model.fit(X, y)

# ----------------------------------------------------------------------
# 3.  Predict probabilities
# ----------------------------------------------------------------------
print("Predicting probabilities …")
proba = model.predict_proba(X)[:, 1]     # P(class = 1)

# ----------------------------------------------------------------------
# 4.  Save augmented master file
# ----------------------------------------------------------------------
X_prob = X.copy()
X_prob["xgb_pos_prob"] = proba

X_prob.to_csv(out_path, index=False)
print(f"Saved → {out_path.resolve()}")
