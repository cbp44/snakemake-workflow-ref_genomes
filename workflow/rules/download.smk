rule download_genome:
    output:
        fasta="resources/{species}/genome.fa",
    conda:
        "../envs/curl.yaml"
    log:
        "logs/download-{species}-genome.log"
    shadow:
        "shallow"
    params:
        release=config['ref']['release'],
    # cache: True
    shell:
        # Note: --insecure used for large files because https fails to connect, but http times out
        "(curl http://ftp.ensembl.org/pub/release-{params.release}/fasta/{wildcards.species}/dna_index/CHECKSUMS | sed 's/  */ /g' | cut -d ' ' -f 3 > files.txt; "
        "curl -o {output.fasta}.gz --insecure https://ftp.ensembl.org/pub/release-{params.release}/fasta/{wildcards.species}/dna_index/$(grep '.fa.gz$' files.txt); "
        "gunzip {output.fasta}.gz) &> {log}"


rule download_annotation:
    output:
        gtf="resources/{species}/genome.gtf",
    conda:
        "../envs/curl.yaml"
    log:
        "logs/download-{species}-annotation.log"
    shadow:
        "shallow"
    params:
        release=config['ref']['release'],
    # cache: True
    shell:
        # Note: --insecure used for large files because https fails to connect, but http times out
        "(curl http://ftp.ensembl.org/pub/release-{params.release}/gtf/{wildcards.species}/CHECKSUMS | sed 's/  */ /g' | cut -d ' ' -f 3 > files.txt; "
        "curl -o {output.gtf}.gz --insecure https://ftp.ensembl.org/pub/release-{params.release}/gtf/{wildcards.species}/$(grep '{params.release}.gtf.gz$' files.txt); "
        "gunzip {output.gtf}.gz) &> {log}"
