import os.path


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
# genomeChrBinNbits == min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])
# genomeSAindexNbases == min(14, log2(GenomeLength)/2 - 1)
rule star_index:
    input:
        fasta="resources/ensembl/genome.fa",
        fai="resources/ensembl/genome.fa.fai",
        annotation="resources/ensembl/genome.gtf",
    output:
        directory("resources/ensembl/star_genome"),
    conda:
        "../envs/star.yaml"
    params:
        sjdbOverhang="99",  # Sequencing read lenegth - 1
        genomeSAindexNbases="12",  # Calculated from genome size, see STAR docs
        outTmpDir=lambda wc, output: os.path.join(output[0], "_STARtmp"),
    log:
        "logs/ensembl/star_index.log",
    threads: 24
    shell:
        "mkdir -p {output}; "
        "STAR --runThreadN {threads} --runMode genomeGenerate"
        " --outTmpDir {params.outTmpDir} --genomeDir {output}"
        " --genomeSAindexNbases {params.genomeSAindexNbases}"
        " --genomeFastaFiles {input.fasta}"
        " --sjdbOverhang {params.sjdbOverhang}"
        " --sjdbGTFfile {input.annotation} {input.annotation} &> {log}"
