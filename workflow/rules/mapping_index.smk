rule genome_faidx:
    input:
        "resources/{species}/genome.fa",
    output:
        "resources/{species}/genome.fa.fai",
    log:
        "logs/{species}-genome_faidx.log",
    wrapper:
        "0.80.2/bio/samtools/faidx"

## TODO: Dynamically adjust genomeSAindexNbases parameter based on size of genome.
rule star_index:
    input:
        fasta="resources/{species}/genome.fa",
        fai="resources/{species}/genome.fa.fai",
        annotation="resources/{species}/genome.gtf",
    output:
        directory("resources/{species}/star_genome"),
    conda:
        "../envs/star.yaml",
    threads: 1
    params:
        extra=lambda wc, input: "--sjdbGTFfile {0} --genomeSAindexNbases 12".format(
            input.annotation)
        ,
        sjdbOverhang="99",  # Sequencing read lenegth - 1
        fasta=lambda wc, input: "--genomeFastaFiles {0}".format(input.fasta),
        gtf=lambda wc, input: "{0}".format(input.annotation),
    log:
        "logs/{species}-star_index.log",
    shell:
        "STAR --runMode genomeGenerate {params.extra} --runThreadN {threads} --genomeDir {output} {params.fasta} --sjdbOverhang {params.sjdbOverhang} {params.gtf} &> {log}"
