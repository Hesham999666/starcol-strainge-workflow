# Donor B StrainGE reference database

This directory contains the workflow used to reproduce the custom
Donor B StrainGE reference database used for strain tracking.

## Reference genomes

The database was constructed from 53 draft bacterial genome assemblies
obtained from cultured Donor B isolates.

Raw whole-genome sequencing reads and draft genome assemblies for these
53 isolates are deposited in the European Nucleotide Archive (ENA):

**Study accession: PRJEB128176**

ENA:
https://www.ebi.ac.uk/ena/browser/view/PRJEB128176

The REFIT2 participant metagenomes are controlled research data and are
not included in this repository.

## Workflow

The reference-database workflow is:

1. Obtain the 53 Donor B draft genome assemblies.
2. Prepare each assembly as a compressed FASTA reference with unique,
   reference-specific contig names.
3. Convert every reference FASTA to a StrainGST HDF5 k-mer set using
   k = 23.
4. Calculate all-vs-all reference similarities using Jaccard, ANI and
   reciprocal subset similarity metrics.
5. Combine all 53 reference k-mer sets into the StrainGST pan-genome
   database.

For 53 references, all-vs-all comparison produces:

    53 × 52 / 2 = 1,378 reference pairs

All 53 isolate references were retained in the final Donor B database.

## Scripts

### `prepare_references.sh`

Prepares draft genome assemblies for use by StrainGE by:

- assigning unique reference-specific contig headers,
- writing sequences in linearized FASTA format,
- compressing the prepared references as `.fa.gz`.

The original project-specific `fasta_adapter` helper script was not
retained. This script is therefore a reconstructed public implementation
of the documented preparation procedure.

### `build_database.sh`

Builds the StrainGST database by:

- running `straingst kmerize -k 23` for each reference,
- calculating all-vs-all similarities with `straingst kmersim`,
- retaining all supplied references,
- creating the final pan-genome database with `straingst createdb`.

The original project-specific database helper script was not retained.
This public script reconstructs the database-building procedure from
the surviving project artifacts and documented StrainGST commands.

## Expected output

For the Donor B reference set:

- 53 prepared `.fa.gz` reference files
- 53 individual `.hdf5` k-mer sets
- `DonorB_similarities.tsv`
- `DonorB_references.txt`
- `DonorB_pan-genome-db.hdf5`

The similarity table contains 1,378 pairwise reference comparisons.
