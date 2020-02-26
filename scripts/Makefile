# Variables
basedir = $(abspath .)

reads = ${basedir}/data/HG001.GRCh38_chr22_0.01_1.fq.gz ${basedir}/data/HG001.GRCh38_chr22_0.01_2.fq.gz
genome = "${basedir}/reference/Homo_sapiens_assembly38_chr22.fa"
annot = "${basedir}/reference/Homo_sapiens_assembly38_chr22.fa"

OUTPUT_DIR = ${basedir}/test-commands/quality-control
fastqc_outputs = ${OUTPUT_DIR}/HG001.GRCh38_chr22_0.01_1_fastqc.html ${OUTPUT_DIR}/HG001.GRCh38_chr22_0.01_1_fastqc.zip ${OUTPUT_DIR}/HG001.GRCh38_chr22_0.01_2_fastqc.html ${OUTPUT_DIR}/HG001.GRCh38_chr22_0.01_2_fastqc.zip

default: show-info fastqc

${OUTPUT_DIR}:
	@mkdir -p $@

.PHONY: show-info
show-info:
	@echo "TEST MAKE - VARIANT ANALYSIS"
	@echo $(abspath .)

${fastqc_outputs} &: ${OUTPUT_DIR} ${reads}
	@fastqc --outdir ${OUTPUT_DIR} ${reads}

.PHONY: fastqc
fastqc: ${fastqc_outputs}
	@echo Finished
