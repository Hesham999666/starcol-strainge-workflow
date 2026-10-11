#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Version 2025-12-17a
// For def. input and output folder see nextflow.config

process readCounter {
    publishDir "${params.outdir}/CallStats", mode: 'copy', overwrite: true

    cpus 1
    memory 4.GB
    time 1.h
    maxForks 19

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_reads_in_fastqs.txt")

    script:
    """
    r1=\$(zcat ${reads[0]} | wc -l)
    r2=\$(zcat ${reads[1]} | wc -l)

    echo \$(( (r1 + r2) / 4 )) > ${sample_id}_reads_in_fastqs.txt
    """
}

process Kmerize {
//    publishDir "${params.outdir}/kmerize", mode: 'copy', overwrite: true
    
    conda "${projectDir}/../environment/workflow.yml"

    cpus 1
    memory 40.GB
    time 1.h
    maxForks 4

    input:
    tuple val(sample_id), path(reads) 

    output:
    tuple val(sample_id), path("${sample_id}.hdf5")

    script:
    """
    straingst kmerize -k 23 -o ${sample_id}.hdf5 ${reads[0]} ${reads[1]}
    """
}

process KmerCompare {
    publishDir "${params.outdir}/KmerCompare", mode: 'copy', overwrite: true

    conda "${projectDir}/../environment/workflow.yml"

    cpus 1
    memory 30.GB
    time 2.h
    maxForks 5

    input:
    tuple val(sample_id), path(hdf5)

    output:
    tuple val(sample_id),
        path("${sample_id}.strains.tsv"),
        path("${sample_id}.stats.tsv"),
        emit: strains

    script:
    """
    straingst run -i ${params.iterations} -s ${params.score} -O -o ${sample_id} ${params.db_dir}${params.db_name}_pan-genome-db.hdf5 ${hdf5}
    """
}

process PrepareRef {
    publishDir "${params.outdir}/PrepareRef", mode: 'copy', overwrite: true

    conda "${projectDir}/../environment/workflow.yml"

    cpus 2
    memory 1.GB
    time 1.h
    maxForks 9

    input:
    tuple val(sample_id), path(strains_tsv), path(_)

    output:
    tuple val(sample_id),
          path("${sample_id}.refs_concat.fasta"),
          path("${sample_id}.refs_concat.meta.json"),
          path("${sample_id}.refs_concat.repetitive.bed"),
          path("${sample_id}.refs_concat.collapsed.tsv"),
          emit: ref_files

    script:
    """
    straingr prepare-ref -s ${strains_tsv} \\
        -p ${params.db_dir}${params.db_name}/{ref}.fa.gz \\
        -S ${params.db_dir}${params.db_name}_similarities.tsv \\
        -o ${sample_id}.refs_concat.fasta
    """
}

process BwaMapping {
//  publishDir "${params.outdir}/BwaMapping", mode: 'copy', overwrite: true

    conda "${projectDir}/../environment/workflow.yml"

    cpus 8
    memory 4.GB
    time 1.h
    maxForks 3

    input:
    tuple val(sample_id),
          path(ref_fasta),
          path(meta_json),
          path(repetitive_bed),
          path(collapsed_tsv),
          path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.bam"),
          path("${sample_id}.bam.bai"),
          emit: bam_pair

    script:
    """
    bwa index ${ref_fasta}

    bwa mem -I 300 -t ${task.cpus} ${ref_fasta} ${reads[0]} ${reads[1]} | samtools sort -@ 2 -O BAM -o ${sample_id}.bam

    samtools index ${sample_id}.bam
    """
}

process CallStats {
    // errorStrategy 'ignore'

    publishDir "${params.outdir}/CallStats", mode: 'copy', overwrite: true

    conda "${projectDir}/../environment/workflow.yml"

    cpus 1
    memory 10.GB
    time 2.d
    maxForks 15

    input:
    tuple val(sample_id), 
        path(ref_fasta), path(meta_json), path(repetitive_bed), path(collapsed_tsv),
        path(bam_file), path(bai_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.called.hdf5"), path("${sample_id}.called.tsv"), emit: callstats_out

    script:
    """
    straingr call ${ref_fasta} ${bam_file} \\
        --hdf5-out ${sample_id}.called.hdf5 \\
        --summary ${sample_id}.called.tsv --tracks all
    """
}


process CallsSummer {
    publishDir "${params.outdir}/CallsSummer", mode: 'copy', overwrite: true
    
    conda "${projectDir}/../environment/workflow.yml"

    cpus 1
    memory 1.GB
    time 1.h
    maxForks 19

    input:
    tuple val(sample_id),
      path(called_hdf5),
      path(called_tsv),
      path(strains_tsv),
      path(stats_tsv),
      path(reads_count_tsv)

    output:
    tuple val(sample_id),
          path("${sample_id}.called_summed.tsv"),
          emit: callsummer_out

    script:
    """
    python3 ${projectDir}/aggregator_calls.py \\
       -strains ${strains_tsv} \\
       -called ${called_tsv} \\
       -reads ${reads_count_tsv} \\
       -o ${sample_id}.called_summed.tsv
    """
}

process CallsSummedPlotter {
    publishDir "${params.outdir}/CallsSummedPlotted", mode: 'copy', overwrite: true
    
    conda "${projectDir}/../environment/workflow.yml"

    cpus 1
    memory 1.GB
    time 1.h
    maxForks 2

    input:
    path(all_summed_tsv)

    output:
    path("CallsSummedPlotted.*")

    script:
    """
    python3 ${projectDir}/plot_calls.py \\
           -callsSummed ${all_summed_tsv} \\
           -o CallsSummedPlotted
    """
}

/* ────────────────────────────────────────────────────────────────
   Workflow definition
──────────────────────────────────────────────────────────────── */

workflow {
      
      paired_reads_ch = Channel.fromFilePairs(
          "${params.input_dir1}/*.fastq.gz",
          { file ->
              def name = file.name

              // Match: anything ... (_|.) (1 or 2) [then optional extra stuff] .fastq.gz
              def m = name =~ /(.+?)[._](1|2)(?:\D.*)?\.fastq\.gz$/

              if (!m) {
                  throw new IllegalArgumentException(
                      "Cannot extract read number (1 or 2) from filename: ${name}"
                  )
            }

            return m[0][1]  // the sample prefix before the final _1/_2/.1/.2
          },
          size: 2
      )

      total_reads     = readCounter(paired_reads_ch)

      kmer_results    = Kmerize(paired_reads_ch)
      strains         = KmerCompare(kmer_results)
      ref_files       = PrepareRef(strains)
      ref_and_reads   = ref_files.join(paired_reads_ch)
      bam_pair        = BwaMapping(ref_and_reads)
      ref_and_bam     = ref_files.join(bam_pair)

      callstats_out   = CallStats(ref_and_bam)
      combinedTables  = callstats_out.join(strains)
      combinedWithReads = combinedTables.join(total_reads)
      
      summedCallsTSV  = CallsSummer(combinedWithReads)

      summed_ch = summedCallsTSV.map { id, tsv -> tsv }

      // Create a single merged TSV file
      merged_tsv_ch = summed_ch.collectFile(
      name: "all_summed.tsv",
      storeDir: "${params.outdir}/CallsSummer",
      keepHeader: true                              // <-- removes duplicate headers!
      )

      CallsSummedPlotter(merged_tsv_ch)
}
