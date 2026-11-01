# Data dictionary and units

The90-row source ledger is a provenance/qualification ledger, not90training samples. observation_id identifies a source-linked paired observation; study_key identifies the recorded paper; source_anchor identifies the original table/figure/page. oat1_uM/oat3_uM and relation fields preserve units and inequalities. endpoint keepsIC50,Ki and apparent-fit families separate. Legacy labels and later overlays are version-attributed and must not be collapsed into a single qualification flag.

The49-row outcomes file has46rows with in_main=True. candidate_structure_id and constitution_smiles are candidate graph keys, not physical reagent certificates. reported_graph_qualified flags the29-row report-associated graph layer. physical_assay_material_certified preserves unresolved status. observed_assay_conditional_delta is log10(IC50_OAT3/IC50_OAT1) after unit alignment. It is not intrinsic affinity, a clinical effect or a binary selective label.

The assay JSON retains three categorical fields and five numeric fields. null means unknown. Probe concentrations useuM, time minutes, temperature Celsius;log transforms and fit-only imputation are defined in the code. Unknown values are not standard laboratory defaults.

The975-prediction table has13procedures evaluated on46and29queries in separate views. Queries and views overlap. predicted_delta and observed_delta are in log10ratio units; absolute_error is their absolute difference. outer_fold is a joint paper-and-constitution component. Equal-query, equal-component and equal-paper summaries target different weighting schemes.

Source-data tables for probe contrasts contain the original endpoint means and uncertainty semantics; later retention tables describe changing query subsets. They do not supply calibrated intervals or independent experimental replicates. Model states and.joblibfiles are saved repository-generated fits; the verifier reconstructs their predictions and compares against retained output records.
