rule Create_Gffutils_DB:
    """Creates a gffutils database that can be used to extract specific feature
    types from the GTF annotation file.
    """
    input:
        "resources/ensembl/genome.gtf.gz"
    output:
        "resources/ensembl/gffutils.db"
    log:
        "logs/create_gffutils_db.log"
    conda:
        "../envs/gffutils.yaml"
    cache: True
    script:
        "../scripts/create_gffutils_db.py"


if is_human_genome():
    use rule Create_Gffutils_DB as Create_Gffutils_DB_MANE with:
        input:
            "resources/ensembl/MANE.GRCh38.v1.0.gtf.gz"
        output:
            "resources/ensembl/mane-gffutils.db"
        log:
            "logs/create_gffutils_db_mane.log"