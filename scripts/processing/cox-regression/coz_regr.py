import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import CoxPHFitter
import numpy as np
import sys

# --- Input --- prelcic cancer file
csv_path = "/Users/tusharsingh/Desktop/motifs_may/may24/output/TRB/shared_motifs/shared_motifs_4k/single/preclic/uterus_trb_survival_ready.csv"
df = pd.read_csv(csv_path)

# Detect motif columns
motif_cols = [col for col in df.columns if col not in ["OS_days", "Status"]]

# Convert motif columns to numeric
for motif in motif_cols:
    df[motif] = pd.to_numeric(df[motif], errors="coerce")

# Cox regression
motif_results = []
for motif in motif_cols:
    sub_df = df[["OS_days", "Status", motif]].dropna()
    if sub_df[motif].sum() < 5 or (sub_df[motif] == 0).sum() < 5:
        continue

    cph = CoxPHFitter()
    try:
        cph.fit(
            sub_df.rename(columns={"OS_days": "duration", "Status": "event"}),
            duration_col="duration",
            event_col="event"
        )
        hr = cph.hazard_ratios_[motif]
        p = cph.summary.loc[motif, "p"]
        motif_results.append({"Motif": motif, "p_value": p, "HR": hr})
    except Exception as e:
        print(f"Skipped {motif}: {e}")

# Prepare results
res_df = pd.DataFrame(motif_results)
res_df["-log10(p)"] = -np.log10(res_df["p_value"])
res_df["Effect"] = res_df["HR"].apply(lambda x: "Protective" if x < 1 else "Risk")
