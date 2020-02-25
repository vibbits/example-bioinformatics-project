#!/usr/bin/env nextflow

//inputData = Channel.create()
//inputPairedReads = Channel.create()
// define threads: 24

// define 1% reads in the data/ folder & reference genome (chr22 limited versions) 
params.reads = "$baseDir/data/HG001.GRCh38_chr22_0.01_{1,2}.fq.gz"
params.genome = "$baseDir/reference/Homo_sapiens_assembly38_chr22.fa"
params.annot = "$baseDir/reference/Homo_sapiens_assembly38.genes_chr22.gff"

/* 
 * prints user convenience 
 */
println "TEST NEXTFLOW - VARIANT ANALYSIS"
println "================================"
println "Pipeline based on https://github.com/nextflow-io/crg-course-nov16/"
println "genome             : ${params.genome}"
println "annotat            : ${params.annot}"
println "reads              : ${params.reads}"

/*
 * get a file object for the given param string
 */
genome_file = file(params.genome)
annotation_file = file(params.annot)


/**
 * Quality control fastq
 */

 process fastqc_reads {
  conda "envs/fastq.yml"
  publishDir OUT_DIR + "fastqc"

  input:
  val data from data_ch

  """
  mkdir test-commands/quality-control
  fastqc --outdir test-commands/quality-control/ data/HG001.GRCh38_chr22_0.01_1.fq.gz data/HG001.GRCh38_chr22_0.01_2.fq.gz
  """
}

/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
read_pairs = Channel.fromFilePairs(params.reads, flat: true)

/*
 * Step 2. Maps a read-pair by using bwa mem 
 */
process mapping {
     
    input:
    file 'genome.index.fa' from genome_file 
    file genome_index from genome_index
    set pair_id, file(read1), file(read2) from read_pairs
 
    output:
    set pair_id, "tophat_out/accepted_hits.bam" into bam
 
    """
    //tophat2 -p ${task.cpus} genome.index ${read1} ${read2}

    mkdir test-commands/bwa-mapping
    bwa mem -t ${task.threads} -M reference/Homo_sapiens_assembly38_chr22.fa data/HG001.GRCh38_chr22_0.01_1.fq.gz data/HG001.GRCh38_chr22_0.01_2.fq.gz | samtools view -b - -o test-commands/bwa-mapping/HG001_chr22_rawmappings.bam
    """
}

/**
 * Quality control bam
 */
 process fastqc_reads {
  conda "envs/fastq.yml"
  publishDir OUT_DIR + "fastqc"

  input:
  val data from data_ch

  """
  fastqc -f bam --outdir test-commands/quality-control/ test-commands/bwa-mapping/HG001_chr22_rawmappings.bam
  """
}

