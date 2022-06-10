rule Download_MANE_Annotation:
    output:
        multiext("resources/ensembl/MANE.GRCh38.v1.0", 
            ".gtf.gz", 
            ".transcripts_by_gene.tsv.gz")
    conda:
        "../envs/curl.yaml"
    log:
        "logs/ensembl/download-mane.log",
    params:
        # --insecure used because https fails cert validation
        # and http times out for larger files
        curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
        vcf=lambda wc: f"https://ftp.ensembl.org/pub/data_files/{config['ref']['species']}/{HUMAN_ACC}/mane/MANE_v1.0/MANE.GRCh38.v1.0.gtf.gz",
        tsv=lambda wc: f"https://ftp.ensembl.org/pub/data_files/{config['ref']['species']}/{HUMAN_ACC}/mane/MANE_v1.0/MANE.GRCh38.v1.0.transcripts_by_gene.tsv.gz",
    cache: True
    retries: 3
    shell:
        """
        (curl -o {output[0]} {params.curl} {params.vcf} \
         && curl -o {output[1]} {params.curl} {params.tsv}) &> {log}
        """