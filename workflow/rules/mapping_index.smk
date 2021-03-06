import os.path


rule Get_Genome_Length:
    """Calculate the total length of the genome from the FASTA sequence.
    """
    input:
        fasta="resources/ensembl/genome.fa.gz",
    output:
        ensure(temp("resources/ensembl/genome.fa.seqlen"), non_empty=True),
    params:
        awk=workflow.source_path("../scripts/sequence_length.awk"),
    conda:
        "../envs/awk.yaml"
    shell:
        "awk -f {params.awk} <(zcat {input.fasta}) > {output}"


rule Calc_STAR_Params:
    input:
        "resources/ensembl/genome.fa.seqlen",
    output:
        ensure(temp("resources/ensembl/genome.star_param.genomeSAindexNbases"), non_empty=True)
    script:
        "../scripts/calc_star_param.py"


rule Mapping_Index:
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
        " --genomeFastaFiles {params.fasta}"
        " --sjdbOverhang {params.sjdbOverhang}"
        " --sjdbGTFfile {params.gtf}; "
        "mv Log.out {output}; "
        "rm -rf {params.fasta} {params.gtf} {params.outTmpDir}"
        ") &> {log}"


if is_human_genome():
    use rule Mapping_Index as Mapping_Index_MANE with:
        input:
            fasta="resources/ensembl/genome.fa.gz",
            annotation="resources/ensembl/MANE.GRCh38.v1.0.gtf.gz",
            genomeSAindexNbases="resources/ensembl/genome.star_param.genomeSAindexNbases",
        output:
            directory("resources/ensembl/star_genome_mane"),
        log:
            "logs/ensembl/star_index_mane.log",