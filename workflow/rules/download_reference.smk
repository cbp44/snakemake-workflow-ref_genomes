rule determine_fasta_file:
    """
    Determine the correct Ensembl genome fasta file to get.

    If there is a dna.primary_assembly.fa.gz file, get that, otherwise get the
    dna.toplevel.fa.gz one.
    """
    output:
        temp("resources/ensembl/genome_fasta_file.txt"),
    conda:
        "../envs/curl.yaml"
    log:
        "logs/ensembl/determine_fasta_file.log",
    params:
        # --insecure used because https fails cert validation
        # and http times out for larger files
        curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
        base_url=lambda wc: get_ensembl_base_url("fasta", "dna"),
        file_ext="-e 'dna.primary_assembly.fa.gz$' -e 'dna.toplevel.fa.gz$'",
    shadow:
        "shallow"
    shell:
        "(curl {params.curl} {params.base_url}CHECKSUMS "
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


rule download_genome_fasta:
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
        base_url=lambda wc: get_ensembl_base_url("fasta", "dna"),
    shadow:
        "shallow"
    shell:
        "("
        "files=$(cat {input}); "
        'curl -o {output} {params.curl} {params.base_url}"$files"'
        ") &> {log}"


rule download_gene_annotation:
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
        base_url=get_ensembl_base_url("gtf"),
        file_ext=lambda wc: "^{0}\..+\.[[:digit:]]+\.gtf\.gz$".format(
            config["ref"]["species"].capitalize()
        ),
    shadow:
        "shallow"
    shell:
        "(curl {params.curl} {params.base_url}/CHECKSUMS "
        "   | sed 's/  */ /g' "
        "   | cut -d ' ' -f 3 > files.txt; "
        "curl -o {output.gtf} {params.curl} {params.base_url}$(grep -E '{params.file_ext}' files.txt); "

        "rm -f files.txt"
        ") &> {log}"
        # "gunzip {output.gtf}.gz; "


# rule download_vcf_annotation:
#     output:
#         vcf="resources/ensembl/{vcf_type}.vcf.gz",
#         csi="resources/ensembl/{vcf_type}.vcf.gz.csi",
#     conda:
#         "../envs/curl.yaml"
#     log:
#         "logs/ensembl/download-vcf-{vcf_type}.log",
#     params:
#         # --insecure used because https fails cert validation
#         # and http times out for larger files
#         curl="--insecure --retry 3 --retry-connrefused --show-error --silent",
#         base_url=lambda wc: get_ensembl_base_url(
#             "variation/vcf",
#             "{0}_{1}".format(wc.species, wc.vcf_type),
#         ),
#     wildcard_constraints:
#         species="homo_sapiens",
#         vcf_type="|".join(
#             [
#                 "clinically_associated",
#                 "phenotype_associated",
#                 "somatic",
#                 "somatic_incl_consequences",
#                 "structural_variations",
#             ]
#         ),
#     shadow:
#         "shallow"
#     shell:
#         "("
#         "   curl -o {output.vcf} {params.curl} {params.base_url}.vcf.gz; "
#         "   curl -o {output.csi} {params.curl} {params.base_url}.vcf.gz.csi"
#         ") &> {log}"
