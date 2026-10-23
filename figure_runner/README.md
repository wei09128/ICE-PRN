# ICE-PRN figures v5.12

Rebuild exactly Figures1–4 and S1–S3 from saved estimates; no model fitting or new statistical tests.

```bash
conda activate ice_prn_cnmf
cd /mnt/c/Wei/ICE-PRN
python ICE_PRN_FIGURES_v5_12/build_all.py
```

Default output: figures/. Optional `--out /path/to/new_folder`. Requires the packages in requirements.txt. The runner validates frozen inputs and writes BUILD_STATUS.json.

Main sequence: transfer, scoring, biological associations, capacity and p53 readouts. Figure3D combines broad cancer-pathway, matched-reference and estrogen contrasts; all signs and inadequate-balance markers are retained. Figure4 promotes the original-model ribosomal contrast; the full secondary grid stays inS3. S1 retains directional p53 coverage; S2 shows absolute program skill without duplicating main relative contrasts. Four main plus three supplementary figures,21 PNG/PDF/SVG files.

Captions match the revised JCIM documents. Current source-data filenames follow current numbering. Old layout and removed paper analyses are archived, not erased from the research record. Original within-cancer-pathway endpoints remain in source_data even though the main panel now displays relative contrasts. Conditional intervals and exploratory endpoint families have not changed. No disease efficacy, pathway mediation, biological replication or direction interaction is established.

WRITING_STYLE.md carries the user's concise-writing preferences. EDITORIAL_CHANGES.md explains the revision. The preceding v5.11 version remains recoverable through file version history and the archived layout records.
