rule CDS_Regions:
    """Pull out the coding sequence regions and save to a tsv file.
    """
    input:
        db="resources/ensembl/gffutils.db",
    output:
        tsv="resources/ensembl/cds_regions.tsv",
    conda:
        "../envs/gffutils.yaml"
    cache: True
    script:
        "../scripts/gtf_to_cds_tsv.py"




if is_human_genome():
    use rule CDS_Regions as CDS_Regions_MANE with:
        input:
            db="resources/ensembl/mane-gffutils.db",
        output:
            tsv="resources/ensembl/mane-cds_regions.tsv",
        log:
            "logs/create_gffutils_db_mane.log"