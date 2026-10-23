# Reproducibility scope

## Included and verified

The v5.12 figure renderer uses saved estimates, intervals and mappings. It rebuilds all seven active figures without large raw-expression matrices. The TOC renderer uses the recorded rounded benchmark skill values. Frozen source/code hashes are in RELEASE_HASHES.json; figure inputs have an additional original INPUT_HASHES.json.

## Included runners with external requirements

Unzip an analysis packet into a separate directory and read its own README and contract before running. Each packet preserves its original directory layout. Matching and readout analyses may require the large LINCS files at user-configured WSL paths. The ribosomal packet contains the returned per-target losses needed to replay that completed follow-up. These packets do not replace an end-to-end raw-data workflow.

The complete genetic/chemical raw-data extraction, original fitting environment and every historical matching run are not reconstructed by the figure renderer. Do not describe a successful figure build as independent replication or as rerunning every analysis from raw data. Large GCTX/H5AD files remain at the original source or on the analysis workstation.

## Release checks

`python verify_release.py` checks source hashes, original frozen figure-input hashes, ZIP integrity, archive paths and the absence of files at or above GitHub's 100 MiB per-file limit. `VALIDATION.json` records the actual successful build and Python/package versions used for this snapshot. No new scientific inference is run during packaging.

Generated PDF/SVG/PNG files, logs and build-status records are excluded from the release source manifest so routine regeneration can change export metadata. Source inputs remain protected by hashes. Run the verifier before rendering.

## Historical material

`figure_runner/archive/` preserves superseded layouts, source indexes and removed analyses. Inputs not drawn in the active seven figures remain for provenance. The old chemical-opposition analysis is not a current manuscript claim. The October 1 v7.26.0 handoff is superseded for this manuscript; it is not installed as an active task list here.
