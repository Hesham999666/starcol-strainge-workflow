# STARCOL StrainGE workflow

Reproducible strain-tracking workflow used for the study:

**Strain-level analysis of donor Bifidobacterium engraftment after faecal microbiota transplantation in irritable bowel syndrome**

This repository contains the computational workflow used to construct a
custom Donor B StrainGE reference database, detect reference-related
strains in shotgun metagenomes, compare recipient strains with the exact
administered donor sample, and classify donor-strain engraftment and
persistence.

## Study overview

The analysis was performed within the REFIT2 faecal microbiota
transplantation (FMT) trial.

A custom StrainGE database was constructed from 53 cultured bacterial
isolates from Donor B.

The complete 53-reference database was used for strain-level profiling.
Downstream analyses in the manuscript focused on four cultured
Bifidobacterium strains:

- Bifidobacterium longum
- Bifidobacterium bifidum
- Bifidobacterium adolescentis
- Bifidobacterium animalis

B. animalis was retained for descriptive analyses but excluded from the
primary engraftment analyses because donor-concordant signals were also
observed in placebo recipients.

## Workflow overview

The computational workflow consists of four main stages:

    Donor B isolate assemblies
              |
              v
    StrainGST reference database
       (53 isolate genomes)
              |
              v
       StrainGE Caller
              |
              v
    strict strain detection
              |
              +----------------------+
              |                      |
              v                      v
       recipient HDF5        exact donor HDF5
              |                      |
              +----------+-----------+
                         |
                         v
                  StrainGE Compare
                         |
                         v
              donor-strain concordance
                         |
                         v
               confirmed engraftment

For T12 samples, the same Caller and exact-donor Compare criteria were
used to assess persistence among strains already confirmed as engrafted
at T03.

---

# 1. Reference genome data

Raw whole-genome sequencing reads and draft genome assemblies for the
53 Donor B bacterial isolates used to construct the StrainGE reference
database are deposited in the European Nucleotide Archive (ENA).

**ENA study accession: PRJEB128176**

https://www.ebi.ac.uk/ena/browser/view/PRJEB128176

REFIT2 participant metagenomic and clinical data are controlled research
data and are not distributed through this repository.

---

# 2. Software

The workflow was developed using:

- StrainGE 1.3.9
- Nextflow 24.04.2
- Python 3.10.19
- BWA 0.7.19
- samtools 1.22.1
- MUMmer4 4.0.1
- SeqKit 2.12.0
- pandas 2.3.3
- NumPy 1.26.4
- matplotlib 3.10.7
- seaborn 0.13.2

The Conda environment is provided in:

    environment/workflow.yml

Create the environment with:

    conda env create -f environment/workflow.yml

and activate it with:

    conda activate starcol-strainge

Nextflow 24.04.2 should also be available on the system.

---

# 3. Building the Donor B reference database

The reference database contains all 53 cultured Donor B isolate genomes.

Prepare the assemblies:

    mkdir -p database/DonorB

    database_build/prepare_references.sh \
        assemblies \
        database/DonorB

The preparation step assigns unique reference-specific contig names,
linearizes the FASTA sequences and compresses them as `.fa.gz`.

Build the StrainGST database:

    cd database

    ../database_build/build_database.sh \
        DonorB \
        DonorB

    cd ..

For each reference, StrainGST k-merization uses:

    k = 23

The 53 references produce:

    53 × 52 / 2 = 1,378

all-vs-all reference comparisons.

Expected database structure:

    database/
    ├── DonorB/
    │   ├── reference_1.fa.gz
    │   ├── reference_1.hdf5
    │   ├── ...
    │   └── reference_53.hdf5
    ├── DonorB_similarities.tsv
    ├── DonorB_references.txt
    └── DonorB_pan-genome-db.hdf5

All 53 references are retained in the final pan-genome database.

See:

    database_build/README.md

for additional details.

---

# 4. StrainGE Caller

Paired-end metagenomic FASTQ files are placed under:

    data/fastq/

The Caller workflow performs:

1. FASTQ read counting
2. StrainGST k-merization
3. StrainGST reference detection
4. StrainGR reference preparation
5. BWA read mapping
6. StrainGR variant calling
7. per-reference aggregation

The core StrainGST detection command uses:

    straingst run \
        -i 50 \
        -s 0.0000002 \
        -O

Run the Caller workflow with:

    nextflow run caller/main_strainGE_caller.nf \
        -profile saga

The combined Caller output is written under:

    results/CallsSummer/all_summed.tsv

---

# 5. StrainGE Compare

StrainGE Compare was used to compare each recipient sample with the
corresponding exact administered donor metagenome.

The pair definition is supplied through a tab-separated table containing:

    comparison_id
    recipient_id
    recipient_timepoint
    donor_id
    donor_hdf5
    recipient_hdf5

A synthetic example is provided in:

    example_metadata/exact_donor_pairs.example.tsv

No real REFIT2 participant identifiers are included in this repository.

Run exact donor-recipient comparisons with:

    compare/run_exact_pairs.sh \
        example_metadata/exact_donor_pairs.example.tsv

The core comparison command is:

    straingr compare donor.called.hdf5 recipient.called.hdf5

---

# 6. Strain classification

## Broad detection

A reference-related signal was considered broadly detected when:

    abundance_metagenome > 0.01%

## Strict Caller detection

Strict strain detection required all of the following:

    abundance_metagenome > 0.01%
    callablePct > 10%
    confirmedPct >= 95%
    snpPct < 1%
    gapCount = 0
    gapLength = 0

## Exact-donor genomic concordance

Donor-recipient genomic concordance required:

    commonPct >= 10%
    singleAgreePct >= 99.95%

## Confirmed T03 engraftment

Confirmed engraftment at T03 required BOTH:

1. strict Caller detection in the recipient T03 sample; and
2. genomic concordance with the exact administered donor sample.

## T12 persistence

Persistence at T12 required:

1. confirmed T03 engraftment;
2. strict Caller detection at T12; and
3. genomic concordance at T12 with the exact administered donor sample.

The classification definitions are provided in:

    classification/criteria.yml

Thresholds can be applied with:

    classification/apply_thresholds.py

Caller and exact-donor Compare results can then be combined using:

    classification/combine_pair.py

and T03-to-T12 persistence can be assessed using:

    classification/call_persistence.py

### Percentage units

The aggregation scripts calculate `Pct` variables as fractions between
0 and 1.

For example:

    10%    = 0.10
    95%    = 0.95
    99.95% = 0.9995

Similarly:

    0.01% metagenomic abundance = 0.0001

The classification scripts use these raw fractional values.

---

# 7. Repository structure

    starcol-strainge-workflow/
    ├── README.md
    ├── nextflow.config
    ├── environment/
    │   └── workflow.yml
    ├── database_build/
    │   ├── README.md
    │   ├── prepare_references.sh
    │   └── build_database.sh
    ├── caller/
    │   ├── main_strainGE_caller.nf
    │   ├── aggregator_calls.py
    │   └── plot_calls.py
    ├── compare/
    │   ├── main_strainGE_compare.nf
    │   ├── aggregator_relations.py
    │   ├── plot_relations.py
    │   └── run_exact_pairs.sh
    ├── classification/
    │   ├── README.md
    │   ├── criteria.yml
    │   ├── apply_thresholds.py
    │   ├── combine_pair.py
    │   └── call_persistence.py
    ├── example_metadata/
    │   └── exact_donor_pairs.example.tsv
    └── docs/

---

# 8. Workflow provenance

The public workflow was prepared from the computational scripts used
during the STARCOL/REFIT2 strain-level analysis.

Project-specific paths, usernames, data-transfer commands and
participant identifiers have been removed from the public versions.

The original project-specific helper scripts used for FASTA preparation
and database construction were not retained. The scripts in
`database_build/` therefore reconstruct those steps from the documented
procedure, surviving database artifacts and StrainGST commands.

The core StrainGE Caller and Compare commands and analysis parameters
are preserved from the surviving workflow.

---

# 9. Data availability

Reference isolate sequencing data:

    ENA PRJEB128176

REFIT2 participant metagenomic and associated clinical data are not
publicly distributed because they contain participant-level research
data. Access is governed by the REFIT2 data-sharing and ethics
requirements.

---

# 10. License

This workflow is released under the MIT License.

See:

    LICENSE

for the full license text.

---

# Citation

A formal citation for the associated manuscript and a Zenodo DOI for
this workflow will be added when available.
