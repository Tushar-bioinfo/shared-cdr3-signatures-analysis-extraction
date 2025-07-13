import pandas as pd
from lifelines import KaplanMeierFitter, statistics
import matplotlib.pyplot as plt
from pathlib import Path


# -------- CONFIG --------
RECEPTOR = "TRB"  # 🔁 Change this to IGH, IGK, etc.

input_dir = Path("/Users/tusharsingh/Desktop/motifs_may/may24/output/TRB/shared_motifs/shared_motifs_4k/single")
output_root = Path(f"/Users/tusharsingh/Desktop/motifs_may/may24/output/TRB/shared_motifs/shared_motifs_4k/single/output")
output_root.mkdir(exist_ok=True)

# -------- GET FILES --------
csv_files = list(input_dir.glob("*_survival_ready.csv"))

for file_path in csv_files:
    cancer_label = file_path.stem.replace("_survival_ready", "")
    print(f"\n📂 Processing: {cancer_label}")

    df = pd.read_csv(file_path)
    motif_cols = [col for col in df.columns if col not in ["Case_ID", "OS_days", "Status"]]

    # Set up output dirs
    cancer_dir = output_root / cancer_label
    plot_dir = cancer_dir / "km_plots"
    cancer_dir.mkdir(exist_ok=True)
    plot_dir.mkdir(exist_ok=True)

    results = []

    for motif in motif_cols:
        try:
            sub_df = df[["OS_days", "Status", motif]].dropna()
            group1 = sub_df[sub_df[motif] == 1]
            group0 = sub_df[sub_df[motif] == 0]

            if len(group1) < 5 or len(group0) < 5:
                with open(cancer_dir / "skipped_motifs.txt", "a") as f:
                    f.write(f"{motif}\tPresent: {len(group1)}\tAbsent: {len(group0)}\n")
                continue  # skip motifs with very small groups

            # Fit KM
            kmf1 = KaplanMeierFitter()
            kmf0 = KaplanMeierFitter()
            kmf1.fit(group1["OS_days"], group1["Status"], label=f"{motif} Present (n={len(group1)})")
            kmf0.fit(group0["OS_days"], group0["Status"], label=f"{motif} Absent (n={len(group0)})")

            # Log-rank test
            result = statistics.logrank_test(group1["OS_days"], group0["OS_days"],
                                             event_observed_A=group1["Status"],
                                             event_observed_B=group0["Status"])
            pval = result.p_value

            results.append({
                "Motif": motif,
                "P_value": pval,
                "n_present": len(group1),
                "n_absent": len(group0)
            })

            if pval < 0.1:
                ax = kmf1.plot(ci_show=True, show_censors=True, ci_alpha=0.2, color="#1B9CFC")
                kmf0.plot(ax=ax, ci_show=True, show_censors=True, ci_alpha=0.2, color="#FC427B")
                plt.title(f"{cancer_label} | {motif} Survival (p={pval:.3f})", fontsize=12)
                plt.xlabel("Days", fontsize=12)
                plt.ylabel("Survival Probability", fontsize=12)
                plt.xticks(size=8)
                plt.yticks(size=8)
                plt.grid(True)
                plt.tight_layout(pad=0.1)
                plt.savefig(plot_dir / f"{motif}_km.png")
                plt.close()

        except Exception as e:
            print(f"⚠️ Skipped {motif} in {cancer_label}: {e}")

    # Save summary
    summary_df = pd.DataFrame(results)
    if not summary_df.empty:
        summary_df.sort_values("P_value", inplace=True)
        summary_df.to_csv(cancer_dir / "motif_survival_summary.csv", index=False)

    print(f"✅ {cancer_label} done — summary + plots saved.")
