# ICE-PRN figures — manuscript v5.5

Rebuilds Figures 1–3 and S1–S4 from bundled saved results. No fitting, downloads or new statistical tests.

```bash
cd /mnt/c/Wei/ICE-PRN/ICE_PRN_FIGURES_v5_5
python -m pip install -r requirements.txt
python build_all.py            # or: python build_all.py --out /path/to/dir
```

| Figure | Content | Code | Was in v4.17 |
|---|---|---|---|
| 1 | Copy/zero vs learned skill, genetic (K562↔RPE1) and chemical (MCF7↔PC3), plus normalized genetic | build_figure1.py | Fig 1A + Table S4 |
| 2 | Comparator and magnitude effects on scores | build_figure2.py | Fig 1B–F |
| 3 | Combined rank/capacity (A), combined p53 coverage/matching (B), identity-program contrasts (C), drug disease-readout coverage (D) | build_figure3.py | S4 + Fig 2A + Table S3 |
| S1 | p53 coverage by direction | build_s1.py | S1 |
| S2 | Identity-TF and drug-mechanism annotations | build_s2.py | Fig 2B–C |
| S3 | Fixed identity-program readout skill | build_s3.py | S3 |
| S4 | Per-chemical MCF7/PC3 reference-opposition patterns and technical-half agreement | build_s4.py | New exploratory analysis |

Removed from this paper (kept in v4.17 package for a separate AD paper): iPSC/MAT2A, Trem2 and D-45113 (old Figs 3, 4, S2).

`_style.py` holds the shared palette (direction blue #1b6fa8 / orange #c55f2c; outcome groups teal, ochre, slate, neutral grey), checked for colour-vision separation and contrast.

`inputs/normalized_cross_gene_skill.tsv` carries the two normalized cross-gene skills and intervals exactly as reported in manuscript v4.16 Methods; the normalization archive itself does not include the result table. Replace it with the archived result file when available.

## Layout update

Figure 2 uses a 2 × 2 arrangement. Panel C contains both directional coverage curves on separate vertical scales; former E is now D. Remove manuscript references to Figure 2E: use Figure 2D. References to the coverage curves are Figure 2C. Overall plot titles and bottom explanatory paragraphs are moved to Figure_captions.md. Short panel headings, axes, keys and numerical annotations are retained. Figure 3 headings describe predictive capacity without claiming biological dimensionality. No numerical inputs, estimates, fits or intervals changed.

## Correction in v5.5

The supplied reference image is Figure 3 (capacity and p53), not Figure 2 (comparator and risk). Its former A+B are now new A, former C+D are now new B; new C shows identity-program contrasts and new D shows drug disease-pathway gene coverage. `build_all.py` rebuilds this corrected composite together with the other five figures. Figure 2 is restored to its original v5.3 builder and layout. The mistaken standalone Figure_2_v6 is superseded and is not produced by this package.

The full original S2 annotation analysis and S3 program results are retained for completeness. No manuscript text or numbering is changed automatically. Manuscript references to former Figure 3 C/D should now point to Figure 3 B; former Figure 3 B (rank-gain share) becomes Figure 3 A. New identity/drug references are Figure 3 C/D. Source tables and the updated caption distinguish program-skill differences in C from gene-annotation coverage differences in D.

## v5.5: exploratory drug figure and complete rebuild

`build_all.py` now produces seven figures (1–3, S1–S4), each as PNG/PDF/SVG: 21 graphics. Figures 1–3 and S1–S3 retain their v5.4 builders, inputs and interpretation. The current Figure 3D disease-pathway annotation endpoint is retained. The new per-chemical opposition analysis is Figure S4, not a replacement drug-efficacy panel.

`build_s4.py` reads frozen per-identity scores from the completed return, reproduces all category counts, and draws both technical halves against two separate prostate references. It does not recompute expression responses, fit models, download external data or run new significance tests. Gene-mask and reference-quality differences limit cross-context score comparisons. No p-values, significance claims, clinical winner, phenotype restoration or therapeutic selectivity are implied by quadrant membership.

Pending before considering a main drug figure: independent viability matching/validation, activity beyond broad toxicity, and a separate held-out test of whether learned transfer predicts destination opposition scores better than copy/mean. These outcomes are not available and are not fabricated as panels. No additional WSL scoring run is required to rebuild the completed figures.
