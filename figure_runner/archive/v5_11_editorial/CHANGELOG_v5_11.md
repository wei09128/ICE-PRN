# v5.11 (23 October 2026)
- Figure 3 now holds biological associations only (A p53 coverage, B p53 matching, C identity programs, D cancer pathways); rank panels moved.
- New Figure 4: A skill by rank, B share of rank-32 gain (both formerly Fig 3A), C p53 target readouts vs other readouts by rank, D p53-minus-other gain over affine (12 contrasts, adjusted 99.583%). Inputs: inputs/rank_p53 from ICE_PRN_RANK_P53_READOUT_v1.
- New Figure S4: estrogen-response readouts (PC3->MCF7) and matched cancer-pathway readouts. Inputs: inputs/readout_followups from ICE_PRN_ESTROGEN_READOUT_TEST_v1 and ICE_PRN_DRUG_READOUT_MATCHING_v1.
- Saved estimates only; no models refitted, no new statistical tests. Rank-curve 'all readouts' values in the rank-p53 return match full_rank_curve.tsv to 1.6e-14.
- Run: python build_all.py
- In-figure explanatory notes removed from Figs 3, 4 and S4; that information is in the legends.
- New Figure S5: ribosomal-target x p53-readout contrast (inputs/ribosome_p53 from ICE_PRN_RIBOSOMAL_P53_TEST_v1).
