
rule Download_VCF_Annotation:
    output:
        "resources/ensembl/{vcf_type}_variants.vcf.gz"
        # multiext("resources/ensembl/{vcf_type}_variants", 
        #     ".vcf.gz", 
        #     ".vcf.gz.csi")
    conda:
        "../envs/curl.yaml"
    log:
        "logs/ensembl/download-vcf-{vcf_type}.log",
    params:
        # --insecure used because https fails cert validation
        # and http times out for larger files
        curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
        base_url=get_ensembl_url("variation/vcf"),
        # vcf_url=lambda wc: get_ensembl_url("variation/vcf", f"{config['ref']['species']}_{wc.vcf_type}.vcf.gz"),
        vcf=lambda wc: f"{config['ref']['species']}_{wc.vcf_type}.vcf.gz",
        # csi_url=lambda wc: f"{config['ref']['species']}_{wc.vcf_type}.vcf.gz.csi",
        # csi_url+lambda wc: get_ensembl_url("variation/vcf", f"{config['ref']['species']}_{wc.vcf_type}.vcf.gz.csi"),
        # server_filename=lambda wc: f"{config['ref']['species']}_{wc.vcf_type}",
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
        # """
        # (curl -o {output[0]} {params.curl} {params.base_url}/{params.vcf} \
        #  && curl -o {output[1]} {params.curl} {params.base_url}/{params.csi_url}) &> {log}
        # """

# TODO: Make this programmatically determine GRCh38 version from genome_info.json
rule Download_ClinVar_Variants:
    # input:
    output:
        "resources/clinvar/clinvar.vcf.gz"
        # multiext("resources/clinvar/clinvar", 
        #          ".vcf.gz",
        #          ".vcf.gz.tbi")
    params:
        curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
        vcf="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    conda: "../envs/curl.yaml"
    log: "logs/ensembl/download-clinvar-vcf.log"
    cache: True
    retries: 3
    shell:
        "curl -o {output[0]} {params.curl} {params.vcf} &> {log}"
        # """
        # (curl -o {output[0]} {params.curl} {params.vcf}; \
        # curl -o {output[0]}.tbi {params.curl} {params.vcf}.tbi) &> {log}
        # """
        


rule Filter_Unusual_Variants:
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