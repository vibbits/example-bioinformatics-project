# Variables
basedir = $(abspath .)

reads = ${basedir}/data/HG001.GRCh38_chr22_0.01_1.fq.gz ${basedir}/data/HG001.GRCh38_chr22_0.01_2.fq.gz
annot = ${basedir}/reference/Homo_sapiens_assembly38_chr22.fa

FASTQC_DIR = ${basedir}/test-commands/quality-control
fastqc_outputs = ${FASTQC_DIR}/HG001.GRCh38_chr22_0.01_1_fastqc.html ${FASTQC_DIR}/HG001.GRCh38_chr22_0.01_1_fastqc.zip ${FASTQC_DIR}/HG001.GRCh38_chr22_0.01_2_fastqc.html ${FASTQC_DIR}/HG001.GRCh38_chr22_0.01_2_fastqc.zip

BWA_DIR = ${basedir}/test-commands/bwa-mapping
bwa_output = ${BWA_DIR}/HG001_chr22_rawmappings.bam

default: show-info all

${FASTQC_DIR} ${BWA_DIR}:
	@mkdir -p $@

.PHONY: show-info
show-info:
	@echo "TEST MAKE - VARIANT ANALYSIS"
	@echo $(abspath .)

${fastqc_outputs} &: ${FASTQC_DIR} ${reads}
	@fastqc --outdir $^

${bwa_output}: ${annot} ${reads}
	@bwa mem -t 4 -M $^ | samtools view -b - -o $@
	@echo Finished bam

.PHONY: all
all: ${FASTQC_DIR} ${BWA_DIR} ${fastqc_outputs} ${bwa_output}
	@fastqc --outdir ${FASTQC_DIR} ${bwa_output}
	@echo Finised workflow
