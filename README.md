# example-bioinformatics-project
Compare and contrast available tools for bioinformatics projects

## Workflow tools
May be based on [this](https://github.com/grst/snakemake_nextflow_wdl).

Some [documentation for CWL](https://www.commonwl.org/).

There are many tools. [An informal survey or popularity](https://docs.google.com/forms/d/e/1FAIpQLScoj8Po4P3Qrh7rbJrq2R35c3PQsNCynEeEVUAdcGyly7TT_Q/viewanalytics).

Some examples of tools we may test are:
* [GNU Make](https://www.gnu.org/software/make/)
* [SnakeMake](https://snakemake.readthedocs.io)
* [Nextflow](https://www.nextflow.io/)
* [Toil](https://toil.ucsc-cgl.org/)
* [Cromwell/WDL](https://cromwell.readthedocs.io)
* [Guix Workflow Language](https://www.guixwl.org/)
* [Apache Taverna]()
* [Apache Airflow](https://airflow.apache.org/)
* [Cuneiform]()
* [Ruffus](http://www.ruffus.org.uk/)

James: GNU Make & Guix Workflow Language
Tuur: Snakemake & Nextflow

Set-up: analyze the data of variant analysis course: quality control & alignment of subset of the data.

## Situation 1
We presume the data & reference have already been downloaded, if not, download it similarly as explained in [this repo](https://github.com/tmuylder/variantcalling/blob/master/Commands_jan2020_10pc.md). The commands that should be used are listed here below.

Respect the following file structure:
```
 project-folder/
    |
    |- data/
    |   |- HG001.GRCh38_chr22_0.01.fq.gz
    |   |- ...
    |
    |- reference/
    |   |- ...
    |
    |- nextflow/ 
    |   |- ...
    |
    |- snakemake/
    |   |-...
    | 
    ...
```

## Setup
A conda environment can be set up with
```bash
conda env create -n workflow -f environment.yml
conda activate workflow
```

## Commands

```
# Make directory quality-controls & perform fastqc
mkdir test-commands/quality-control
fastqc --outdir test-commands/quality-control/ data/HG001.GRCh38_chr22_0.01_1.fq.gz data/HG001.GRCh38_chr22_0.01_2.fq.gz

# Make directory bwa-mappings & perform alignment & convert sam to bam
mkdir test-commands/bwa-mapping
bwa mem -t 4 -M reference/Homo_sapiens_assembly38_chr22.fa data/HG001.GRCh38_chr22_0.01_1.fq.gz data/HG001.GRCh38_chr22_0.01_2.fq.gz | samtools view -b - -o test-commands/bwa-mapping/HG001_chr22_rawmappings.bam

# Quality control on bam file
fastqc --outdir test-commands/quality-control/ test-commands/bwa-mapping/HG001_chr22_rawmappings.bam
```
## Evaluation and Experiences
### GUIX Workflow Language
This tool is very new. Based on GUIX which is very promising for reproducability (minimal system bootstrap). The language is based on high-level definitions of tasks and workflows. It is a DSL for scheme so a user could drop into a powerful general-purpose programming language if they need to.

However I could not even install it. It just had a `0.2.0` release 1 week ago and the mailing list has bugs reported that will be fixed in a `0.2.1` release "soon". In summary it is not ready to use _yet_.
