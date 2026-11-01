# OAT1 and OAT3 selectivity evidence and assay context

Version1.0.0 contains the data and code supporting an exploratory evaluation of paired OAT1/OAT3 inhibition selectivity under chemical distribution shift. It is distinct from the author's published single-target OAT1 docking benchmark (doi:10.1016/j.insi.2026.100707); the old docking dataset DOI must not be used to identify this release.

## Scientific scope

The source ledger has90rows. The mixed analysis has46ordinary pairedIC50records from12papers and43candidate constitutions; the report-associated graph layer has29records from8papers. These overlapping layers are not additive. Graph association does not certify physical reagent identity. Ki, censored, apparent-fit and unknown states are preserved. No independent external validation or interval calibration was performed.

Mixed-layer chemical-family selection MAE was0.572801 versus training mean0.623067. Retraining within the29-record layer gave0.715729 versus mean0.654749. Assay increments depended on the view, weighting, missingness and scaling controls. All families and unfavorable cases are retained. These are exploratory findings, not a deployed or externally validated predictor.

## Reproduction

Use Python3.12 with the pinned packages. From this repository root run:

```sh
python -m pip install -r requirements.txt
python scripts/check_identity_assay_increment.py
```

The numerical replay checks all975predictions, exact source/constitution exclusions, fit-only encoders, parameter-selection rules, linear/kernel decompositions, manual forest traversal and matched comparisons without fitting new models. Saved.joblibfiles are repository-generated forest states. The ROOT marker and historical-package-only guard changes are documented in PORTABILITY_CHANGES.csv.

To rerun the fixed13-procedure nested comparison on the same inputs, make a working copy of this repository and run:

```sh
python scripts/train_identity_assay_increment.py
python scripts/check_identity_assay_increment.py
python scripts/plot_identity_assay_increment.py
```

The training script writes regenerated outputs into the included run directory. Do not edit the fixed grids or interpret a rerun as an independent test. Reproducing the exact original forests requires the recorded scikit-learn version and seed. The smaller output checks can run withoutRDKit; descriptor generation in the earlier adequacy script usesRDKit.

## Files and qualification

The runs directory preserves the relative scientific input layout: source ledger, mixed-identity outcomes with inclusion flags, raw assay metadata, ten descriptor values, Tanimoto matrix, nested partitions, inner scores, outer states and all predictions. INPUT_SHA256.csv records the public input bytes, while RELEASE_SHA256.csv records release files. Metadata paths referring to private source archives are generalized; published DOI and table/page anchors remain. Source PDFs, publisher full text, credentials and host-specific directories are excluded.

The former110-observation/27-source response inventory is not a complete footprint of historical abstract visibility. Source qualification and reused literature should be assessed separately from reproducibility of a saved prediction. This release makes no claim of pristine prospective blinding or original wet-lab measurement.

## Rights and citation

Original code is released underMIT. Author-created curation annotations and derived outputs useCC BY4.0; pre-existing source measurements, citations and database identifiers retain their original attribution and rights. This license does not relicense third-party articles or figures. SeeTHIRD_PARTY_NOTICES.md and CITATION.cff. This repository contains data/code, not the unpublished manuscript or the author's supplied published PDF.
