# Example Bioinformatics Project
Compare and contrast available tools for bioinformatics projects

## Workflow tools
May be based on [this](https://github.com/grst/snakemake_nextflow_wdl).

Some [documentation for CWL](https://www.commonwl.org/).

There are many tools. [An informal survey or popularity](https://docs.google.com/forms/d/e/1FAIpQLScoj8Po4P3Qrh7rbJrq2R35c3PQsNCynEeEVUAdcGyly7TT_Q/viewanalytics).

Some examples of tools we may test are:
* [GNU Make](https://www.gnu.org/software/make/) :heavy_check_mark:
* [SnakeMake](https://snakemake.readthedocs.io) :construction:
* [Nextflow](https://www.nextflow.io/) :heavy_check_mark:
* [Toil](https://toil.ucsc-cgl.org/)
* [Cromwell/WDL](https://cromwell.readthedocs.io)
* [Guix Workflow Language](https://www.guixwl.org/) :x:
* [Apache Taverna]()
* [Apache Airflow](https://airflow.apache.org/) :heavy_check_mark:
* [Cuneiform]()
* [Ruffus](http://www.ruffus.org.uk/)


Results summary:

:heavy_check_mark: :heavy_check_mark: - recommended

:heavy_check_mark: - usable

:x: - not (yet) useful

:construction: - evaluation in progress

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
This tool is very new. Based on GUIX which is very promising for reproducability (minimal system bootstrap). The language is based on high-level definitions of "processes" and "workflows". It is a DSL for scheme so a user could drop into a powerful general-purpose programming language if they need to.

The simple sample workflow from the documentation looks like this:
```
process greet
  packages "hello"
  # { hello }

process sleep
  packages "coreutils"
  # {
    echo "Sleeping..."
    sleep 10
  }

process eat (with something)
  name
    string-append "eat-" something
  # {
    echo "Eating {{something}}"
  }

process bye
  # { echo "Farewell, world!" }

workflow simple-wisp
  processes
    let
      :
        eat-fruit : eat "fruit"
        eat-veges : eat "vegetables"
      graph
        eat-fruit -> greet
        eat-veges -> greet
        sleep     -> eat-fruit eat-veges
        bye       -> sleep
```

I could install gwl, but the documented `guix workflow` command does not exist after installation. It just had a `0.2.1` bugfix release which does not appear to fix this issue. In summary, it is not ready to use _yet_. The [manual](https://www.guixwl.org/manual/gwl.html) seems fairly comprehensive but this has, of course, not been tried in anger.

### GNU Make
This tool is _very_ old and also very common on GNU/Linux systems (users may not even need to install it). The workflow is expressed as a tree of dependency relationships between files. Independent trees may be executed in parallel. Make will not re-run tasks that _do not need to be re-run_ where this is defined by the timestamp on the output relative to the input. These relationships are defined like this:
```make
target: input
 recipe
multiple targets: several dependencies
 recipe $^ > $@ # Generate a single target
many outputs &:: that are all
 generated at once by this recipe
.PHONY: just-a-name
just-a-name:
 recipe
```
[Documentation](https://www.gnu.org/software/make/manual/html_node/index.html) is good and comprehensive though lacking in examples. Though this is generally made up for because `make` is used so foten. A quick internet search will reveal many potential answers, unfortunately some of which are well out of date. That being said, `make` doesn't change very much. There are 2 main downsides to expressing workflows in `make` (besides the modern features that other, more special purpose tools have) are the arcane syntax and the difficult-to-reason-about semantics of dependency trees described manually. For software, `Makefile`s are usually automatically generated.

### Nextflow
Nextflow is on the way to be the reference bioinformatics pipeline tool. The scripting language is Groovy which makes it in my opinion daunting at first, however with some time spent reading the [manual](https://www.nextflow.io/docs/latest/getstarted.html) it becomes more comprehensible. Together with a very troublesome error-reporting, these are Nextflow's main drawbacks. It has many strengths on the other side: 
- Easy installation and running
- Use of conda en containers for environment with all necessary tools, so not caring about any tools that need to be installed (+ reproducibility). 
- Pipeline resumes after it breaks down at a specific point
- Theoretically: parallelization,  scalability

Script so far is [here](scripts/testmapping.nf)

IMO: I think there is a lack on specific bioinformatics related tutorials although it is used specifically for this purpose. E.g.: according to the manual, `set` is deprecated in the current version (hence no more documentation on how to use it) and has been replaced by `tuple`, however all of the pipelines out there use `set`. `tuples` should be generated from channels with read pairs, however it's continuously getting stuck somewhere and it's close to impossible to find out where the error is exactly. 

### Snakemake
TBC.
Might be interesting to check: 
- https://github.com/maxplanck-ie/snakepipes/tree/2.0.0
- https://vincebuffalo.com/blog/2020/03/04/understanding-snakemake.html
- https://github.com/vanheeringen-lab/snakemake-workflows

### Apache Airflow
Initial setup is very complex. Pipelines are just Python scripts. Has a visualisation of the pipeline as it runs. Very general purpose, leading to the extra complexity. Documentation seems ok, but frastratingly not quite what one is looking for in many cases. For example, the tutorial states, "Please take the time to understand **<concept>**" but does not provide _any_ pointers to documentation for that concept. Feels verbose and "enterprisy".
 ![Enterprise software](https://i.imgflip.com/138d0a.jpg)
