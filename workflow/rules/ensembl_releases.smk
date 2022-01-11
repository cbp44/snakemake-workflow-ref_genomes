

rule get_ensembl_genome_list:
    output:
        "results/ensembl/release-{release}/available_genomes.txt",
    log:
        "logs/ensembl/release-{release}/get_genome_list.txt",
    conda:
        "../envs/curl.yaml"
    shadow:
        "shallow"
    shell:
        "(curl -o {output} http://ftp.ensembl.org/pub/release-{wildcards.release}/species_EnsemblVertebrates.txt) &> {log}"


rule extract_supported_species:
    input:
        "results/ensembl/release-{release}/available_genomes.txt",
    output:
        "results/ensembl/release-{release}/supported_genomes.txt",
    params:
        supported_species="|".join(config.get("supported_species", ["mus_musculus"])),
    conda:
        "../envs/awk.yaml"
    log:
        "logs/ensembl/release-{release}/extract_supported_species.txt",
    shadow:
        "shallow"
    shell:
        "((head -n 1 {input} | sed 's/^#//'; grep -E {params.supported_species:q} {input}) | cut -f  > {output}) &> {log}"
