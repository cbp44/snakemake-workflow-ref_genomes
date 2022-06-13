rule Transcript_BED:
    """Creates a BED file containing all of the annotated transcript sites
    for use by RESeQC in the RNA-seq workflow to determine coverage of reads.
    """
    input:
        db="resources/ensembl/gffutils.db",
    output:
        bed="resources/ensembl/transcript_annotation.bed",
    log:
        "logs/rseqc_gtf2bed.log"
    conda:
        "../envs/gffutils.yaml"
    cache: True
    script:
        "../scripts/gtf_to_transcript_bed.py"
