rule Download_Genome_Metadata:
    output:
        json="resources/ensembl/{metadata_type}.json",
    params:
        extra=lambda wc, output: "--species {0} {2} --output {1}".format(config["ref"]["species"], output.json, wc.metadata_type),
    conda:
        "../envs/rest.yaml"
    wildcard_constraints:
        metadata_type="|".join(["genome_info","variation_consequences"])
    shadow:
        "shallow"
    cache: True
    retries: 3
    script:
        "../scripts/ensembl_rest_api.py"


rule Determine_FASTA_File:
    """
    Determine the correct Ensembl genome fasta file to get.

    If there is a dna.primary_assembly.fa.gz file, get that, otherwise get the
    dna.toplevel.fa.gz one.
    """
    output:
        ensure(temp("resources/ensembl/genome_fasta_file.txt"), non_empty=True),
    conda:
        "../envs/curl.yaml"
    log:
        "logs/ensembl/determine_fasta_file.log",
    params:
        # --insecure used because https fails cert validation
        # and http times out for larger files
        curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
        base_url=lambda wc: get_ensembl_url("fasta", "dna"),
        file_ext="-e 'dna.primary_assembly.fa.gz$' -e 'dna.toplevel.fa.gz$'",
    shadow:
        "shallow"
    cache: True
    retries: 3
    shell:
        "(curl {params.curl} {params.base_url}/CHECKSUMS "
        "   | sed 's/  */ /g' "
        "   | cut -d ' ' -f 3 | grep {params.file_ext} > files.txt; "
        "cat files.txt; "
        "set +e; grep primary_assembly files.txt; "
        "if [ $? -eq 0 ]; "
        "then grep -e 'dna.primary_assembly.fa.gz$' files.txt > {output} ; "
        "else grep -e 'dna.toplevel.fa.gz$' files.txt > {output}; "
        "fi; set -e; "
        "rm -f files.txt; "
        ") &> {log}"


rule Download_Genome_FASTA:
    input:
        "resources/ensembl/genome_fasta_file.txt",
    output:
        "resources/ensembl/genome.fa.gz",
    conda:
        "../envs/curl.yaml"
    log:
        "logs/ensembl/download-genome.log",
    params:
        # --insecure used because https fails cert validation
        # and http times out for larger files
        curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
        base_url=lambda wc: get_ensembl_url("fasta", "dna"),
    shadow:
        "shallow"
    cache: True
    retries: 3
    shell:
        "("
        "files=$(cat {input}); "
        'curl -o {output} {params.curl} {params.base_url}/"$files"'
        ") &> {log}"


rule Download_Gene_Annotation:
    output:
        gtf="resources/ensembl/genome.gtf.gz",
    conda:
        "../envs/curl.yaml"
    log:
        "logs/ensembl/download-annotation.log",
    params:
        # --insecure used because https fails cert validation
        # and http times out for larger files
        curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
        base_url=get_ensembl_url("gtf"),
        file_ext=lambda wc: "^{0}\..+\.[[:digit:]]+\.gtf\.gz$".format(
            config["ref"]["species"].capitalize()
        ),
    shadow:
        "shallow"
    cache: True
    retries: 3
    shell:
        "(curl {params.curl} {params.base_url}/CHECKSUMS "
        "   | sed 's/  */ /g' "
        "   | cut -d ' ' -f 3 > files.txt; "
        "curl -o {output.gtf} {params.curl} {params.base_url}/$(grep -E '{params.file_ext}' files.txt); "
        "rm -f files.txt"
        ") &> {log}"


rule Download_VCF_Annotation:
    output:
        multiext("resources/ensembl/{vcf_type}_variants", 
            ".vcf.gz", 
            ".vcf.gz.csi")
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
        vcf_url=lambda wc: f"{config['ref']['species']}_{wc.vcf_type}.vcf.gz",
        csi_url=lambda wc: f"{config['ref']['species']}_{wc.vcf_type}.vcf.gz.csi",
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
        """
        (curl -o {output[0]} {params.curl} {params.base_url}/{params.vcf_url} \
         && curl -o {output[1]} {params.curl} {params.base_url}/{params.csi_url}) &> {log}
        """
