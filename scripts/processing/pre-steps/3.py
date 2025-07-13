import pandas as pd
from pathlib import Path

# -------- CONFIG --------
RECEPTOR = "IGH"  # 🔁 Change to "TRB", "IGK", etc. as needed

motif_matrix_dir = Path("/Users/tusharsingh/Desktop/motifs_may/may24/output/IGH/shared_motifs/4mer/motif_matrix")
clinical_dir = Path("/Users/tusharsingh/Desktop/motifs_may/may24/output/TRB/clinical")
output_dir = Path("/Users/tusharsingh/Desktop/motifs_may/may24/output/IGH/shared_motifs/4mer/preclic")
output_dir.mkdir(exist_ok=True)

# Map of motif matrix file prefix to matching clinical file
clinical_files = {
    "kidney_igh": "kidney_clinical.tsv",
    "uterus_igh": "uterus_clinical.tsv",
    "lung_ac_igh": "lung_ac_clinical.tsv",
    "lung_scc_igh": "lung_scc_clinical.tsv",
}

for prefix, clin_file in clinical_files.items():
    print(f"\n📂 Merging for: {prefix}")

    motif_file_name = f"{prefix}_{RECEPTOR}_shared_top20_motifs.csv"
    motif_path = motif_matrix_dir / motif_file_name
    motif_df = pd.read_csv(motif_path)

    dups_motif = motif_df.duplicated(subset="Case_ID").sum()
    motif_df = motif_df.drop_duplicates(subset="Case_ID")
    print(f"🧼 Removed {dups_motif} duplicate Case_IDs from motif matrix")

    # Load clinical data
    clin_df = pd.read_csv(clinical_dir / clin_file, sep="\t")

    # Standardize case ID
    clin_df["Case_ID"] = clin_df["case_submitter_id"].astype(str)

    # Convert survival columns to numeric
    clin_df["days_to_death"] = pd.to_numeric(clin_df["days_to_death"], errors="coerce")
    clin_df["days_to_last_follow_up"] = pd.to_numeric(clin_df["days_to_last_follow_up"], errors="coerce")

    # Compute OS_days and status
    clin_df["OS_days"] = clin_df["days_to_death"].fillna(clin_df["days_to_last_follow_up"])
    clin_df["Status"] = clin_df["vital_status"].map({"Dead": 1, "Alive": 0})

    # Keep only necessary columns
    survival_df = clin_df[["Case_ID", "OS_days", "Status"]].dropna()
    dups_clin = survival_df.duplicated(subset="Case_ID").sum()
    survival_df = survival_df.drop_duplicates(subset="Case_ID")
    print(f"🧼 Removed {dups_clin} duplicate Case_IDs from clinical file")

    # Merge with motif matrix
    merged_df = pd.merge(motif_df, survival_df, on="Case_ID", how="inner")

    # Final deduplication check
    dups_final = merged_df.duplicated(subset="Case_ID").sum()
    if dups_final > 0:
        print(f"⚠️ Still found {dups_final} duplicate Case_IDs in merged result — removing.")
        merged_df = merged_df.drop_duplicates(subset="Case_ID")

    # Save
    out_file = output_dir / f"{prefix}_survival_ready.csv"
    merged_df.to_csv(out_file, index=False)
    print(f"✅ Saved: {out_file} with {merged_df.shape[0]} unique patients")
