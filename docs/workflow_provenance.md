# Workflow provenance

This repository contains a public, cleaned version of the computational
workflow accompanying the STARCOL/REFIT2 StrainGE strain-level analysis.

## Components derived from surviving project code

The following components were derived from workflow files used during
the analysis:

- StrainGE Caller Nextflow workflow
- StrainGE Compare Nextflow workflow
- Caller aggregation script
- Compare aggregation script
- Caller and Compare plotting scripts
- exact donor-recipient pairwise comparison logic
- StrainGE parameters and software versions

Project-specific absolute paths, usernames, data-transfer commands,
storage locations and participant identifiers were removed from the
public versions.

Small changes were made to improve portability and clarity, including:

- replacing absolute paths with repository-relative paths;
- using a single documented Conda environment;
- removing internal data-transfer logic;
- removing the duplicated `.summary.tsv` filename suffix in the public
  Compare workflow;
- retaining the explicit `ref` field in the public Compare summary.

These changes do not alter the strain-detection or genomic-concordance
thresholds used in the study.

## Reconstructed components

The original project-specific helper scripts used for reference FASTA
preparation and StrainGST database construction were not retained.

The following scripts are therefore reconstructed public implementations:

- `database_build/prepare_references.sh`
- `database_build/build_database.sh`

They reproduce the documented database-building procedure using the
surviving database artifacts, documented commands and StrainGST tools.

The Donor B reference database used in the study:

- contained 53 cultured isolate references;
- used StrainGST k-mer size k = 23;
- contained all 1,378 possible pairwise comparisons among the 53
  references;
- retained all 53 references in the final pan-genome database.

## Exact donor mapping

For the final donor-recipient genomic comparison analysis, each recipient
sample was compared with the exact administered Donor B metagenome.

The participant-specific pair table used in the study is not distributed
because it contains controlled study information.

The public repository instead provides:

`example_metadata/exact_donor_pairs.example.tsv`

using synthetic identifiers and example file paths.

## Classification

Confirmed T03 donor-strain engraftment required both:

1. strict StrainGE Caller detection in the recipient T03 sample; and
2. StrainGE Compare genomic concordance with the exact administered
   donor sample.

T12 persistence required confirmed T03 engraftment followed by strict
Caller detection and exact-donor genomic concordance at T12.

The thresholds are documented in:

`classification/criteria.yml`
