# Integration of NMR mQTLs with Cardiovascular Endpoints

This section provides codes for the methods section `Integration of metabolomic
measurements with cardiovascular endpoints`. It consists of a co-localization
analysis of shared genetic variants between NMR metabolites and ASCVD endpoints
(`locus` effects) and a 2-Sample Mendelian Randomization Analysis of the effects
of variants modulating metabolite concentrations on ASCVD outcomes (`level`
effects).

## Colocalization Analysis

The co-localization analysis was conducted using a `nextflow` pipeline that can
be found at
[`https://github.com/comp-med/nf-nmr-coloc`](https://github.com/comp-med/nf-nmr-coloc).
It is added as a submodule under `02_nmr_ascvd_colocalization/`.

## Mendelian Randomization Analysis

All scripts necessary to run the analysis can be found in the directory
`01_nmr_ascvd_mendelian_randomization`. The R script is launched with SLURM and
runs for each of the 249 metabolites.
