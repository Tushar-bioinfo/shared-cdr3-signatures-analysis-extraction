import pandas as pd
from pathlib import Path
from collections import defaultdict
import json
# -------- CONFIG --------
motif_csv_path = "/Users/tusharsingh/Desktop/motifs_may/may24/output/IGH/IGH_k4/Top50_Motifs_Observation_Summary_k4.csv"
cutoffs = [5, 10, 15, 20, 25]
output_dir = Path("/Users/tusharsingh/Desktop/motifs_may/may24/output/IGH/shared_motifs/4mer/top_shared")
output_dir.mkdir(exist_ok=True)

# -------- LOAD --------
df = pd.read_csv(motif_csv_path)

# Get all cancer labels from the column names (skip spacer cols)
cancer_labels = [col for col in df.columns if not col.startswith("Spacer_")]

# -------- Extract Shared Motifs --------
shared_motif_sets = {}
ranked_motifs = defaultdict(list)

for label in cancer_labels:
    ranked = df[label].dropna().tolist()
    ranked_motifs[label] = ranked  # full top 50 per cancer

for k in cutoffs:
    print(f"\n🧩 Shared motifs in Top {k}:")
    top_k_sets = [set(ranked_motifs[label][:k]) for label in cancer_labels]
    shared = set.intersection(*top_k_sets)
    shared_motif_sets[f"shared_top{k}"] = sorted(shared)

    if shared:
        for m in sorted(shared):
            print(f"  ✓ {m}")
    else:
        print("  (none shared in top {})".format(k))

# -------- Save Shared Motifs --------
for k, motif_list in shared_motif_sets.items():
    with open(output_dir / f"{k}_motifs.txt", "w") as f:
        for motif in motif_list:
            f.write(motif + "\n")


with open(output_dir / "shared_motif_sets.json", "w") as f:
    json.dump(shared_motif_sets, f, indent=2)

print("\n✅ Done. All shared motif sets saved to:", output_dir)
