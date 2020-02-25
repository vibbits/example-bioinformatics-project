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
 home/
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
