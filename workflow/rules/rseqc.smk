rule RSeQC_Create_DB:
    input:
        "resources/ensembl/genome.gtf.gz"
    output:
        "results/qc/rseqc/annotation.db"
    log:
        "results/logs/rseqc_create_db.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/create_rseqc_db.py"

# rule RSeQC_GTF2Bed:
#     input:
#         get_gtf_annotation_file()
#     output:
#         bed="results/qc/rseqc/annotation.bed",
#         db="results/qc/rseqc/annotation.db",
#     log:
#         "results/logs/rseqc_gtf2bed.log"
#     conda:
#         "../envs/gffutils.yaml"
#     script:
#         "../scripts/gtf_to_transcript_bed.py"