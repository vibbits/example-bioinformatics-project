#!/usr/bin/env nextflow

/**
 * Parameters & (home) directories
 */
params.abspath = "/home/tuur/example-bioinformatics-project"
params.genome = "$params.abspath/reference/Homo_sapiens_assembly38_chr22.fa" //genome_path
params.outdir = "$params.abspath/nextflow"
params.reads = "$params.abspath/data/*0.01_{1,2}.fq" // reads_ch
params.env = "$params.abspath/environment.yml"

/*
 * prints user convenience
 */
println """\
TEST NEXTFLOW - VARIANT ANALYSIS
================================
Pipeline based on https://github.com/nextflow-io/crg-course-nov16/ & https://github.com/nextflow-io/gccbosc18-training
Reads            : $params.reads
Output-folder    : $params.outdir
Reference        : $params.genome
Environment      : $params.env
"""

/**
 * Set channels
 */
 Channel
         .fromFilePairs(params.reads) // standard option for paired-end reads
         .ifEmpty { 'No read pairs found' } // if no reads: raise error
         .into { read_pairs1_ch; read_pairs2_ch } // set two channels: fastqc and mapping

genome_file = file(params.genome)

/**
 * Quality control fastq
 */
process fastqc_raw_reads {
  conda "$params.env"
  publishDir "$params.outdir/quality-control", mode: 'copy', overwrite: true

  input:
  set sample_id, file(reads) from read_pairs1_ch

  script:
  """
  mkdir -p $params.outdir/quality-control
  fastqc --outdir $params.outdir/quality-control ${reads}
  """
}

/*
 * Step 2. Maps a read-pair by using bwa mem
 */
process mapping {
    conda "$params.env"
    publishDir "${params.outdir}/bwa-mapping", mode: 'copy', overwrite: true

    input:
    set sample_id, file(reads) from read_pairs2_ch

    output:
    set sample_id, file("${sample_id}_raw_mappings.bam") into bam_ch

    script:
    // CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    // readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    """
    bwa mem -M ${genome_file} ${reads[0]} ${reads[1]} | samtools view -b - -o ${sample_id}_raw_mappings.bam
    """
}

/**
 * Quality control bam
 */
process fastqc_maps {
  conda "$params.env"
  publishDir "$params.outdir/quality-control-mapped", mode: 'copy', overwrite: true

  input:
  set sample_id, file(maps) from bam_ch

  """
  mkdir -p $params.outdir/quality-control-mapped
  fastqc -f bam --outdir $params.outdir/quality-control-mapped ${maps}
  """
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

