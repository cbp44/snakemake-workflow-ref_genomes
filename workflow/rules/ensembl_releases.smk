rule get_ensembl_genome_list:
    output:
        "resources/ensembl/release-{release}/available_genomes.txt",
    log:
        "logs/ensembl-{release}-supported_genomes.txt",
    conda:
        "../envs/curl.yaml",
    shell:
        "(curl -o {output} http://ftp.ensembl.org/pub/release-{wildcards.release}/species_EnsemblVertebrates.txt) &> {log}"

rule extract_supported_genomes:
    input:
        "resources/ensembl/release-{release}/available_genomes.txt",
    output:
        "resources/ensembl/release-{release}/supported_genomes.txt",
    params:
        supported_genomes="|".join(config.get('supported_organisms', ['mus_musculus'])),
    conda:
        "../envs/awk.yaml",
    log:
        "logs/ensembl-{release}-supported_genomes.txt",
    shell:
        "((head -n 1 {input} | sed 's/^#//'; grep -E {params.supported_genomes:q} {input}) | cut -f  > {output}) &> {log}"
