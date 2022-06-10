rule Create_Gffutils_DB:
    input:
        "resources/ensembl/genome.gtf.gz"
    output:
        "resources/ensembl/gffutils.db"
    log:
        "logs/rseqc_create_db.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/create_gffutils_db.py"

rule RSeQC_GTF2Bed:
    input:
        db="resources/ensembl/gffutils.db",
    output:
        bed="resources/ensembl/transcript_annotation.bed",
    log:
        "logs/rseqc_gtf2bed.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf_to_transcript_bed.py"
