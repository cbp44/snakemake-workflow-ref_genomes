import os.path
from math import log2


rule get_genome_length:
    input:
        fasta="resources/ensembl/genome.fa.gz",
    output:
        "resources/ensembl/genome.fa.seqlen",
    params:
        awk=workflow.source_path("../scripts/seqlen.awk"),
    conda:
        "../envs/awk.yaml"
    shell:
        "awk -f {params.awk} <(zcat {input.fasta}) > {output}"


rule calc_star_param:
    input:
        "resources/ensembl/genome.fa.seqlen",
    output:
        ensure("resources/ensembl/genome.star_param.genomeSAindexNbases", non_empty=True)
    script:
        "../scripts/calc_star_param.py"


rule star_index:
    input:
        fasta="resources/ensembl/genome.fa.gz",
        annotation="resources/ensembl/genome.gtf.gz",
        genomeSAindexNbases="resources/ensembl/genome.star_param.genomeSAindexNbases",
    output:
        directory("resources/ensembl/star_genome"),
    conda:
        "../envs/star.yaml"
    log:
        "logs/ensembl/star_index.log",
    params:
        sjdbOverhang="99",  # Sequencing read lenegth - 1
        genomeSAindexNbases=lambda wc, input: read_file_line(input.genomeSAindexNbases),
        outTmpDir=lambda wc, output: os.path.join(output[0], "_STARtmp"),
        gtf=lambda wc, input: input.annotation.replace(".gtf.gz", ".gtf"),
        fasta=lambda wc, input: input.fasta.replace(".fa.gz", ".fa"),
    threads: 24
    cache: True
    shell:
        "(mkdir -p {output}; "
        "gunzip -f --keep {input.fasta}; "
        "gunzip -f --keep {input.annotation}; "
        "STAR --runThreadN {threads} --runMode genomeGenerate"
        " --outTmpDir {params.outTmpDir} --genomeDir {output}"
        " --genomeSAindexNbases {params.genomeSAindexNbases}"
        # " --genomeSAindexNbases $(cat {input.genomeSAindexNbases})"
        " --genomeFastaFiles {params.fasta}"
        " --sjdbOverhang {params.sjdbOverhang}"
        " --sjdbGTFfile {params.gtf}; "
        "rm -rf {params.fasta} {params.gtf} {params.outTmpDir}"
        ") &> {log}"
