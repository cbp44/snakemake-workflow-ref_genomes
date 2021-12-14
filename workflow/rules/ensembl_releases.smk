rule get_ensembl_genome_list:
    output:
        "resources/ensembl/release-{release}/available_genomes.txt"
    conda:
        "../envs/curl.yaml"
    shell:
        "curl -o {output[0]} http://ftp.ensembl.org/pub/release-{wildcards.release}/species_EnsemblVertebrates.txt"

rule extract_supported_genomes:
    input:
        "resources/ensembl/release-{release}/available_genomes.txt"
    output:
        "resources/ensembl/release-{release}/supported_genomes.txt"
    params:
        supported_genomes="|".join(config.get('supported_organisms', ['mus_musculus']))
    shell:
        "(head -n 1 {input[0]} | sed 's/^#//'; grep -E {params.supported_genomes:q} {input[0]}) | cut -f  > {output[0]}"
