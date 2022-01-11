import os.path
from math import log2


rule get_genome_length:
    input:
        fasta="resources/ensembl/genome.fa.gz",
    output:
        temp("resources/ensembl/genome.fa.seqlen"),
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
        temp("resources/ensembl/genome.star_param.genomeSAindexNbases"),
    run:
        with open(input[0]) as fin:
            with open(output[0], "w") as fout:
                fout.write(
                    "{}".format(round(min(14, (log2(int(fin.readline())) / 2) - 1)))
                )


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
        outTmpDir=lambda wc, output: os.path.join(output[0], "_STARtmp"),
        gtf=lambda wc, input: input.annotation.replace(".gtf.gz", ".gtf"),
        fasta=lambda wc, input: input.fasta.replace(".fa.gz", ".fa"),
    threads: 24
    shell:
        "(mkdir -p {output}; "
        "gunzip --keep {input.fasta}; "
        "gunzip --keep {input.annotation}; "
        "STAR --runThreadN {threads} --runMode genomeGenerate"
        " --outTmpDir {params.outTmpDir} --genomeDir {output}"
        " --genomeSAindexNbases $(cat {input.genomeSAindexNbases})"
        " --genomeFastaFiles {params.fasta}"
        " --sjdbOverhang {params.sjdbOverhang}"
        " --sjdbGTFfile {params.gtf}; "
        "rm -rf {params.fasta} {params.gtf} {params.outTmpDir}"
        ") &> {log}"
