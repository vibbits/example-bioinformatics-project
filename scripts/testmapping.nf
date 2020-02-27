#!/usr/bin/env nextflow

/**
 * Parameters & (home) directories
 */
params.abspath = "/home/tuur/example-bioinformatics-project"
params.outdir = "${params.abspath}/nextflow"
params.threads = 4


/*
 * Define channels
 */
// define 1% reads in the data/ folder & reference genome (chr22 limited versions)
reads_ch = Channel.fromPath( "${params.abspath}/data/HG001.GRCh38_chr22_0.01_{1,2}.fq.gz" )
genome_path = "${params.abspath}/reference/Homo_sapiens_assembly38_chr22.fa"
Channel.fromFilePairs( "${params.abspath}/data/HG001.GRCh38_chr22_0.01_{1,2}.fq.gz").set{ read_pairs }


/*
 * prints user convenience
 */
println "TEST NEXTFLOW - VARIANT ANALYSIS"
println "================================"
println "Pipeline based on https://github.com/nextflow-io/crg-course-nov16/"
println "Absolute path                   : ${params.abspath}"
println "Output-folder                   : ${params.outdir}"
println "Analyzing the following datasets: ${reads_ch}"
println "Reference                       : ${genome_path}}"



/**
 * Quality control fastq
 */

 process fastqc_raw_reads {
  conda "${params.abspath}/environment.yml"
  publishDir "${params.outdir}/quality-control", mode: 'copy', overwrite: true

  input:
  file raw_reads from reads_ch

//  output:
//  file "*.{zip,html,txt}" into fastqc_results

  """
  mkdir -p ${params.outdir} && mkdir -p ${params.outdir}/quality-control
  fastqc --outdir ${params.outdir}/quality-control/ ${raw_reads}
  """
}


/*
 * Step 2. Maps a read-pair by using bwa mem
 */
process mapping {
    conda "${params.abspath}/environment.yml"
    publishDir "${params.outdir}/bwa-mapping", mode: 'copy', overwrite: true
    
    input:
    file reference from genome_path
    tuple val(x), file('read1'), file('read2') from read_pairs


//    output:
//    set pair_id, "tophat_out/accepted_hits.bam" into bam

    """
    echo Processing mapping reads ${x}
    mkdir -p ${params.outdir}/bwa-mapping
    bwa mem -t ${task.threads} -M ${reference} read1 read2 | samtools view -b - -o ${params.outdir}/bwa-mapping/HG001_chr22_rawmappings.bam
    """
}
