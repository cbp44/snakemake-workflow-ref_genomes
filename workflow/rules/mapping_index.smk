rule genome_faidx:
    input:
        "resources/ensembl/genome.fa",
    output:
        "resources/ensembl/genome.fa.fai",
    log:
        "logs/ensembl/genome_faidx.log",
    wrapper:
        "0.80.2/bio/samtools/faidx"


## TODO: Dynamically adjust genomeSAindexNbases parameter based on size of genome.
rule star_index:
    input:
        fasta="resources/ensembl/genome.fa",
        fai="resources/ensembl/genome.fa.fai",
        annotation="resources/ensembl/genome.gtf",
    output:
        directory("resources/ensembl/star_genome"),
    conda:
        "../envs/star.yaml"
    threads: 1
    params:
        extra=lambda wc, input: "--sjdbGTFfile {0} --genomeSAindexNbases 12".format(
            input.annotation
        ),
        sjdbOverhang="99",  # Sequencing read lenegth - 1
        fasta=lambda wc, input: "--genomeFastaFiles {0}".format(input.fasta),
        gtf=lambda wc, input: "{0}".format(input.annotation),
    log:
        "logs/ensembl/star_index.log",
    shell:
        "STAR --runMode genomeGenerate {params.extra} --runThreadN {threads} --genomeDir {output} {params.fasta} --sjdbOverhang {params.sjdbOverhang} {params.gtf} &> {log}"
