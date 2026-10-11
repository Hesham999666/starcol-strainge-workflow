#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Version 2025-12-03a
// params.input = "/home/tkl001/nextflow_testing/results/CallStats/*.hdf5"

process compareHdf5 {
//    errorStrategy 'ignore'

    publishDir "${params.outdir}/Hdf5Compared", mode: 'copy', overwrite: true

    conda "${projectDir}/../environment/workflow.yml"

    input:
        tuple path(h1), path(h2)

    output:
        path "*.summary.tsv", emit: compareHdf5_out

    script:
    def id1 = h1.getSimpleName().replace('.hdf5','')
    def id2 = h2.getSimpleName().replace('.hdf5','')
    def out = "${id1}_vs_${id2}"

    """
    straingr compare $h1 $h2 -o "${out}.summary.tsv"
    """
}

process ComparisonSummer {
    publishDir "${params.outdir}/ComparisonSummed", mode: 'copy', overwrite: true

    conda "${projectDir}/../environment/workflow.yml"

    input:
        path(summary_tsv)

    output:
        path "*.called_summed.tsv"

    script:
    // derive output filename based on input file name
    def base = summary_tsv.simpleName.replace('.summary.tsv','')
    def out = "${base}.called_summed.tsv"

    """
    python3 ${projectDir}/aggregator_relations.py \
        -compared "${summary_tsv}" \
        -o "${out}"
    """
}

process ConcatSummaries {

    publishDir "${params.outdir}/ComparisonSummed", mode: 'copy', overwrite: true

    input:
        path summed_files

    output:
        path "all_summed.tsv"

    script:
    def first = summed_files[0]
    def rest  = summed_files.size() > 1 ? summed_files[1..-1] : []

    """
    head -n 1 ${first} > all_summed.tsv

    for f in ${rest.join(' ')}; do
        tail -n +2 "\$f" >> all_summed.tsv
    done
    """
}

process ComparisonPlotter {
    publishDir "${params.outdir}/ComparisonPlots", mode: 'copy', overwrite: true

    conda "${projectDir}/../environment/workflow.yml"

    input:
    path(all_summed_tsv)

    output:
    path("ComparisonPlotted.*")

    script:
    """
    python3 ${projectDir}/plot_relations.py \\
           -compSummed ${all_summed_tsv} \\
           -o ComparisonPlotted
    """
}

workflow {

    Channel
        .fromPath("${params.input_dir2}/*.hdf5")
        .collect()
        .flatMap { files ->
            def pairs = []
            for (int i = 0; i < files.size(); i++) {
                for (int j = i+1; j < files.size(); j++) {
                    pairs << [ files[i], files[j] ]
                }
            }
            return pairs
        }
        .set { pairs_ch }

    // step 1: run compare
    summaries_ch = compareHdf5(pairs_ch)

    // step 2: run summer
    summed_ch = ComparisonSummer(summaries_ch)

    // Create a single merged TSV file
    merged_tsv_ch = summed_ch.collectFile(
    name: "all_summed.tsv",
    storeDir: "${params.outdir}/ComparisonSummed",
    keepHeader: true                              // <-- removes duplicate headers!
    )

    // step 3: concatenate all collected files
    //collected_ch = ConcatSummaries(summed_ch.collect())

    // step 4: plot compared data
    ComparisonPlotter(merged_tsv_ch)
}
