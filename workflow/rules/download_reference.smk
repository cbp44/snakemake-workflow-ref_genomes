rule download_genome_fasta:
    output:
        fasta="resources/{species}/genome.fa",
    conda:
        "../envs/curl.yaml"
    log:
        "logs/download-{species}-genome.log",
    shadow:
        "shallow"
    params:
        # --insecure used because https fails cert validation
        # and http times out for larger files
        curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
        base_url=lambda wc: get_ensembl_base_url(
                config["ref"]["release"],
                "fasta",
                wc.species,
                "dna_index"
            ),
        file_ext=".fa.gz",
    cache: True
    shell:
        "(curl {params.curl} {params.base_url}/CHECKSUMS "
        "   | sed 's/  */ /g' "
        "   | cut -d ' ' -f 3 > files.txt; "
        "curl -o {output.fasta}.gz {params.curl} {params.base_url}/" "$(grep '{params.file_ext}$' files.txt)" "; "
        "gunzip {output.fasta}.gz"
        ") &> {log}"


rule download_gene_annotation:
    output:
        gtf="resources/{species}/genome.gtf",
    conda:
        "../envs/curl.yaml"
    log:
        "logs/download-{species}-annotation.log",
    shadow:
        "shallow"
    params:
        # --insecure used because https fails cert validation
        # and http times out for larger files
        curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
        base_url=lambda wc: get_ensembl_base_url(
            config["ref"]["release"],
            "gtf",
            wc.species
        ),
        file_ext="{0}.gtf.gz".format(config["ref"]["release"]),
    cache: True
    shell:
        "(curl {params.curl} {params.base_url}/CHECKSUMS "
        "   | sed 's/  */ /g' "
        "   | cut -d ' ' -f 3 > files.txt; "
        "curl -o {output.gtf}.gz {params.curl} {params.base_url}/" "$(grep '{params.file_ext}$' files.txt)" "; "
        "gunzip {output.gtf}.gz"
        ") &> {log}"


rule download_vcf_annotation:
    output:
        vcf="resources/{species}/{vcf_type}.vcf.gz",
        csi="resources/{species}/{vcf_type}.vcf.gz.csi",
    conda:
        "../envs/curl.yaml"
    log:
        "logs/download-{species}-vcf-{vcf_type}.log",
    shadow:
        "shallow"
    params:
        # --insecure used because https fails cert validation
        # and http times out for larger files
        curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
        base_url=lambda wc: get_ensembl_base_url(
            config["ref"]["release"],
            "variation/vcf",
            wc.species,
            "{0}_{1}".format(wc.species, wc.vcf_type)
        ),
    cache: True
    wildcard_constraints:
        species="homo_sapiens",
        vcf_type="|".join(
            [
                "clinically_associated",
                "phenotype_associated",
                "somatic",
                "somatic_incl_consequences",
                "structural_variations",
            ]
        ),
    shell:
        "("
        "   curl -o {output.vcf} {params.curl} {params.base_url}.vcf.gz; "
        "   curl -o {output.csi} {params.curl} {params.base_url}.vcf.gz.csi"
        ") &> {log}"
