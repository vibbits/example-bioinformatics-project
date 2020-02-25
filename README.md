# example-bioinformatics-project
Compare and contrast available tools for bioinformatics projects

## Workflow tools
May be based on [this](https://github.com/grst/snakemake_nextflow_wdl).
Some [documentation for CWL](https://www.commonwl.org/).
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


## Commands

```
# Make directory quality-controls & perform fastqc
fastqc data/HG001.GRCh38_chr22_0.01.fq.gz data/HG001.GRCh38_chr22_0.01.fq.gz 

# Make directory bwa-mappings & perform alignment & convert sam to bam
bwa mem -t 24 -M -R '@RG\tID:HG001\tLB:NA12878_giab\tPU:unknown-0.0\tPL:Illumina\tSM:NA12878' reference/Homo_sapiens_assembly38_chr22.fa data/HG001.GRCh38_chr22_0.1_1.fq.gz data/HG001.GRCh38_chr22_0.1_2.fq.gz | samtools view -b - -o bwa_mappings/HG001_chr22_rawmappings.bam
```
