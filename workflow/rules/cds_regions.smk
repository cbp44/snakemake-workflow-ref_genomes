rule Genome_CDS_Regions:
    input:
        db="resources/ensembl/gffutils.db",
    output:
        tsv="resources/ensembl/cds_regions.tsv",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf_to_cds_tsv.py"
