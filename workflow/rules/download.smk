# rule get_annotation:
#     output:
#         "resources/genome.gtf.gz",
#     log:
#         "logs/get_annotation.log"
#     params:
#         species=config["ref"]["species"].lower(),
#         species_capitalized=config["ref"]["species"].capitalize(),
#         fmt="gtf",
#         build=config["ref"]["build"],
#         release=config["ref"]["release"]
#     cache: True
#     conda:
#         "../envs/curl.yaml"
#     shell:
#         "curl -o {output[0]} ftp://ftp.ensembl.org/pub/release-{params.release}/gtf/{params.species}/{params.species_capitalized}.{params.build}.{params.release}.gtf.gz 2> {log}"

# rule get_genome:
#     output:
#         "resources/genome.fa.gz",
#     log:
#         "logs/get-genome.log",
#     params:
#         species=config["ref"]["species"],
#         species_capitalized=config["ref"]["species"].capitalize(),
#         datatype="dna",
#         filetype="primary_assembly",
#         build=config["ref"]["build"],
#         release=config["ref"]["release"],
#     cache: True
#     conda:
#         "../envs/curl.yaml"
#     shell:
#         "curl -o {output[0]} ftp://ftp.ensembl.org/pub/release-{params.release}/fasta/{params.species}/{params.datatype}/{params.species_capitalized}.{params.build}.{params.datatype}.{params.filetype}.fa.gz 2> {log}"

rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "0.77.0/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        "logs/get_annotation.log",
    wrapper:
        "0.77.0/bio/reference/ensembl-annotation"
