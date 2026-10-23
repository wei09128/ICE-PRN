# Current figure source map

This index follows the v5.12 layout. `Figure_captions.md` defines the endpoints and intervals. Exact saved tables and mapping records remain in `inputs/` and `source_data/`; frozen input hashes are checked by `build_all.py`.

| Figure | Source inputs | Endpoint |
|---|---|---|
| 1 | benchmark_CI.tsv; inputs/drug; normalized_cross_gene_skill.tsv | Genetic and chemical held-out transfer skill |
| 2 | risk_CI.tsv; uncertainty_summary.json; saved benchmark values | Comparator and response-scale dependence |
| 3A–B | coverage intervals; p53_sensitivity_evidence.zip; recovery intervals | External A549 p53 annotations and matching |
| 3C | inputs/identity_program | Fixed external readout-program skill contrasts |
| 3D | inputs/drug_pathway_skill; readout_followups | Cancer-pathway and estrogen readout contrasts |
| 4 | inputs/rank; rank_p53; ribosome_p53 | Capacity and p53-readout gain, with ribosomal-target follow-up |
| S1 | Saved directional coverage counts and intervals | Directional p53 coverage |
| S2 | inputs/identity_program | Absolute program skill |
| S3 | inputs/rank_p53; ribosome_p53 | Secondary full rank grids |

The manuscript and SI give full bibliographic citations. Genetic responses derive from Replogle et al.; A549 p53 annotations from Lü et al.; chemical responses from LINCS GSE92742; readout memberships from the archived external program/KEGG/Fischer p53 sources. All original source snapshots and mapping records included with the analysis packets retain their original attribution.

The former source index is preserved at archive/SOURCES_predecessor.md. Its figure numbers and chemical-opposition description are historical and do not describe the current paper.
