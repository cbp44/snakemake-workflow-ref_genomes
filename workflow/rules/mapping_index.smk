rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "0.80.2/bio/samtools/faidx"


## TODO: Dynamically adjust genomeSAindexNbases parameter based on size of genome.
rule star_index:
    input:
        fasta="resources/genome.fasta",
        fai="resources/genome.fasta.fai",
        annotation="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: 4
    params:
        extra="--sjdbGTFfile resources/genome.gtf --sjdbOverhang 100ls --genomeSAindexNbases 12",
    log:
        "logs/star_index_genome.log",
    cache: True
    wrapper:
        "0.80.2/bio/star/index"
