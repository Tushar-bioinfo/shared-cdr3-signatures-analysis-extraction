import pandas as pd
import json
from pathlib import Path
from collections import defaultdict

# -------- CONFIG --------
RECEPTOR = "IGH"  # 🔁 Change this to the receptor you are analyzing: "TRB", "IGK", etc.

trb_dir = Path("/Users/tusharsingh/Desktop/raw_data/IGH/run/IGH_motif")
shared_motif_path = Path("/Users/tusharsingh/Desktop/motifs_may/may24/output/IGH/shared_motifs/4mer/top_shared/shared_motif_sets.json")
output_dir = Path("/Users/tusharsingh/Desktop/motifs_may/may24/output/IGH/shared_motifs/4mer/motif_matrix")
output_dir.mkdir(exist_ok=True)

# Load shared motif sets
with open(shared_motif_path) as f:
    shared_motif_sets = json.load(f)

# ✅ Collect all unique motifs from shared_top5 to shared_top20
all_shared_motifs = set()
for k in [5, 10, 15, 20]:
    key = f"shared_top{k}"
    all_shared_motifs.update(shared_motif_sets.get(key, []))
all_shared_motifs = sorted(all_shared_motifs)

print(f"✅ Total motifs collected across top5–20: {len(all_shared_motifs)}")

# -------- Process Each Cancer Type --------
trb_files = list(trb_dir.glob("*.csv"))

for file_path in trb_files:
    print(f"\n📂 Processing {file_path.name}")
    df = pd.read_csv(file_path)
    df = df[df["receptor"] == RECEPTOR]
    df["Case ID"] = df["Case ID"].astype(str).str.split(",").str[0]

    # Group CDR3s per patient
    cdr3s_per_patient = defaultdict(list)
    for _, row in df.iterrows():
        cid = row["Case ID"]
        cdr3 = str(row["CDR3"])
        if isinstance(cdr3, str):
            cdr3s_per_patient[cid].append(cdr3)

    # Build patient x motif matrix
    records = []
    for cid, cdr3s in cdr3s_per_patient.items():
        joined = " ".join(cdr3s)
        record = {"Case_ID": cid}
        for motif in all_shared_motifs:
            record[motif] = int(motif in joined)
        records.append(record)

    df_out = pd.DataFrame(records)
    out_file = output_dir / f"{file_path.stem}_{RECEPTOR}_shared_top20_motifs.csv"
    df_out.to_csv(out_file, index=False)
    print(f"✅ Saved: {out_file}")
