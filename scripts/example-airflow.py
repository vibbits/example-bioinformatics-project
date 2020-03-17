from airflow import DAG
from airflow.operators.bash_operator import BashOperator
from airflow.models import Variable
from airflow.utils.dates import days_ago

from datetime import timedelta
from pathlib import Path


base_path = Variable.get("example_bio_project_base")
out_path = Path(base_path) / 'airflow'
genome_file = Path(base_path) / 'reference' / 'Homo_sapiens_assembly38_chr22.fa'
reads = [Path(base_path) / 'data' / 'HG001.GRCh38_chr22_0.01_1.fq.gz',
         Path(base_path) / 'data' / 'HG001.GRCh38_chr22_0.01_2.fq.gz']
bwa_output = out_path / 'bwa-mapping' / 'HG001_chr22_rawmappings.bam'

default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'start_date': days_ago(0),
    'email': ['james.collier@vib.be'],
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 1,
    'retry_delay': timedelta(minutes=5),
}

dag = DAG(
    'example_bio_workflow',
    default_args=default_args,
    description='Example bioinformatics project',
    schedule_interval=None,
)


raw_reads_command = """
mkdir {{ params.out_path }}/quality-control
fastqc --outdir {{ params.out_path }}/quality-control/ {{ params.reads | join(' ') }}
"""
rawreads = BashOperator(
    task_id='fastq_raw_reads',
    bash_command=raw_reads_command,
    params={'reads': reads, 'out_path': out_path},
    dag=dag,
)

mapping_command = """
mkdir {{ params.out_path }}/bwa-mapping
bwa mem -M {{ params.genome_file }} {{ params.reads | join(' ') }} | samtools view -b - -o {{ params.bwa_output }}
"""
mapping = BashOperator(
    task_id='mapping',
    bash_command=mapping_command,
    params={'genome_file': genome_file, 'reads': reads, 'bwa_output': bwa_output, 'out_path': out_path},
    dag=dag,
)

fastqc_maps_command = """
mkdir {{ params.out_path }}/quality-control-mapped
fastqc -f bam --outdir {{ params.out_path }}/quality-control-mapped {{ params.bwa_output }}
"""
fastqc_maps = BashOperator(
    task_id='fastqc_maps',
    bash_command=fastqc_maps_command,
    params={'out_path': out_path, 'bwa_output': bwa_output},
    dag=dag,
)

[rawreads, mapping] >> fastqc_maps
