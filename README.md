# CDR3 Motif Survival

A Python-based pipeline for identifying shared high-frequency CDR3 amino acid (AA) motifs from tumor-infiltrating lymphocyte repertoires and evaluating their association with overall survival across multiple cancer types.

This project was developed to support the analysis presented in [Pan-cancer, high frequency TRB CDR3 amino acid motifs are associated with outcomes specifically for uterine cancer<img width="468" height="78" alt="image" src="https://github.com/user-attachments/assets/109b5571-791d-4b4c-8335-e069db508dd8" />
], focused on short motif recurrence in TRB and IGH repertoires from CPTAC-3 tumor RNA-seq datasets.

---

## Overview

This repository contains all scripts and workflows required to:
- Extract k-mer motifs (4–11mers) from tumor-derived CDR3 AA sequences.
- Identify shared motifs across cancer types based on top-frequency ranks.
- Construct patient-level binary motif matrices.
- Integrate clinical survival data (OS time, vital status).
- Perform Kaplan–Meier and Cox proportional hazards regression analyses.
- Visualize motif enrichment using heatmaps and survival curves.

---

## Directory Structure 
(coming soon)

