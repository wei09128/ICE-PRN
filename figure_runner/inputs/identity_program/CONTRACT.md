# Frozen external-program readout evaluation

20 October 2026. Exploratory evaluation on already exposed K562/RPE1 data. Program-level prediction outcomes have not been inspected in preparing this package. The fixed program family and mapping masks precede their evaluation.

## Annotation and coverage

Exactly four previously chosen MSigDB snapshots are retained: WELCH_GATA1_TARGETS (10 measured readouts), HALLMARK_HEME_METABOLISM (99), HALLMARK_APICAL_JUNCTION (74), HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION (39). Exact saved symbols only. The local legacy mapping has rows for all 7,226 readouts; two distinct Ensembl IDs map to HSPA14. These two coordinates are excluded from program/complement comparisons, retained in the all-readout protected benchmark. No alias expansion or outcome-guided mapping correction. A current annotation-database release is not established by this export; all conclusions are conditional on the saved mapping. Disagreements with the saved target-symbol table are reported, not silently overwritten.

GATA1-responsive is a mouse-source program with curated human mappings, not a comprehensive human regulator list. Heme is not uniquely erythroid, junction is not RPE1-specific, and EMT is a transition/state program. Do not combine these as an identity score or infer identity causality. Ten GATA1 readouts provide a narrow selected endpoint; no minimum-size exclusion is chosen after seeing outcomes.

## Protected replay and scoring

Reuse unchanged fit/cross_fit helpers and the protected target-grouped fivefold assignments. Four methods: destination mean, source copy, diagonal affine and diagonal-plus-rank-16 cross-gene model. Two directions, two training-half orientations; heldout predictions scored against the opposite destination half. No new model family, hyperparameter selection or model selection. Check the two original response archive hashes and all 80 fold/method construct-loss comparisons against the original reference before releasing aggregate program scores; failure blocks interpretation.

For each program and its complement, MSE is averaged across that mask's readouts excluding the perturbed target's own coordinate when present, then across constructs within target, equally across the 1,874 targets and across the two technical orientations. Complement means all other unambiguously mapped measured readouts, not genes known to lack program activity. Program overlaps are retained and reported. The all-readout protected score includes all 7,226 genes and must reproduce 12.0407% and 9.9424% learned skill within numerical tolerance.

Report model MSE, destination-mean MSE and skill = 1 − model MSE / destination-mean MSE together. Native-control abundance and technical repeatability distributions are descriptive covariates, not an expression/reliability adjustment. Program-versus-complement contrasts are not matched for abundance or measurement precision and cannot attribute performance differences to identity.

## Contrasts and uncertainty

The primary exploratory family comprises exactly eight learned-model contrasts: program skill minus complement skill, four programs × two directions. Preserve both positive and negative results. Ten thousand target-cluster bootstrap draws (seed 2026102010), using the same weights across methods, masks and directions. Report ordinary 95% intervals and approximate Bonferroni 99.375% marginal intervals for these eight contrasts. Other method scores/95% intervals are descriptive; no additional confirmatory tests or P values. No direction heterogeneity test or program-specific gene ranking is authorized by this contract.

Intervals condition on masks, predictions, controls, technical split and prior exploratory choices; they omit retraining, rematching, shared-control uncertainty and independent biological replication. Readout genes are not bootstrap replication units. The target-specific mask losses are exported for audit. Zero/invalid denominators block rather than selectively omit a contrast. No program is enlarged after outcome inspection.

Candidate discovery, CARNIVAL, further p53 exploration and Jurkat outcomes remain parked. No manuscript or figure changes are implied by merely running this exploratory analysis.
