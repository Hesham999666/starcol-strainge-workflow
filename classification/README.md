# Strain classification

The study used two complementary StrainGE steps.

## 1. StrainGE Caller

A reference-related signal was considered broadly detected when:

`abundance_metagenome > 0.01%`

Strict strain detection required all of the following:

- abundance_metagenome > 0.01%
- callablePct > 10%
- confirmedPct >= 95%
- snpPct < 1%
- gapCount = 0
- gapLength = 0

## 2. StrainGE Compare

Genomic concordance with the exact administered donor sample required:

- commonPct >= 10%
- singleAgreePct >= 99.95%

## Confirmed T03 engraftment

Confirmed T03 engraftment required both:

1. strict Caller detection in the recipient T03 sample; and
2. genomic concordance with the exact administered donor sample.

## T12 persistence

Persistence required confirmed T03 engraftment followed by:

1. strict Caller detection at T12; and
2. genomic concordance at T12 with the exact administered donor sample.

## Percentage units

The aggregation scripts in this repository calculate variables such as
`callablePct`, `confirmedPct`, `snpPct`, `commonPct`, and
`singleAgreePct` as fractions between 0 and 1.

For example:

- manuscript callablePct > 10% corresponds to raw value > 0.10
- manuscript confirmedPct >= 95% corresponds to raw value >= 0.95
- manuscript singleAgreePct >= 99.95% corresponds to raw value >= 0.9995

`abundance_metagenome` is also stored as a fraction, so 0.01% corresponds
to a raw value of 0.0001.
