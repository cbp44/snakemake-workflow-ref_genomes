
rule Index_VCF:
    """Creates an index for a given VCF file.
    """
    input:
        "resources/{vcf_file}.vcf.gz"
    output:
        "resources/{vcf_file}.vcf.gz.csi"
    conda:
        "../envs/bcftools.yaml"
    cache: True
    threads: 2
    shell:
        "bcftools index --threads {threads} -o {output} {input}"

if is_human_genome():        
    rule Fetch_VCF_Annotation:
        """Downloads an Ensembl VCF annotation file
        """
        output:
            "resources/ensembl/{vcf_type}_variants.vcf.gz"
        conda:
            "../envs/curl.yaml"
        log:
            "logs/ensembl/download-vcf-{vcf_type}.log",
        params:
            # --insecure used because https fails cert validation
            # and http times out for larger files
            curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
            base_url=get_ensembl_url("variation/vcf"),
            vcf=lambda wc: f"{config['ref']['species']}_{wc.vcf_type}.vcf.gz",
        wildcard_constraints:
            vcf_type="|".join(
                [
                    "clinically_associated",
                    "phenotype_associated",
                ]
            ),
        cache: True
        retries: 3
        shell:
            "curl -o {output[0]} {params.curl} {params.base_url}/{params.vcf} &> {log}"

    rule Fetch_ClinVar_Variants:
        output:
            "resources/clinvar/clinvar.vcf.gz"
        params:
            curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
            vcf=f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{HUMAN_ACC}/clinvar.vcf.gz"
        conda: "../envs/curl.yaml"
        log: "logs/ensembl/download-clinvar-vcf.log"
        cache: True
        retries: 3
        shell:
            "curl -o {output[0]} {params.curl} {params.vcf} &> {log}"
        

    rule Filter_Variants:
        """Remove unusual variant types that cause errors in bedtools.
        """
        input:
            "resources/{vcf_file}.vcf.gz",
            "resources/{vcf_file}.vcf.gz.csi",
        output:
            "resources/{vcf_file}.filtered.vcf.gz"
        conda:
            "../envs/bcftools.yaml"
        cache: True
        threads: 2
        shell:
            "bcftools view -O z9 -o {output} --threads {threads} -V other {input}"
