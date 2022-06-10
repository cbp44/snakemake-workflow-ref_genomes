HUMAN_ACC = "GRCh38"

def get_ensembl_url(filetype, *args):
    """Get a full Ensembl download URL for the configured species given a
    specific filetype (like fasta, gtf, etc.) and optional postfix strings.
    """
    species = config["ref"]["species"]
    base_url = f"ftp://ftp.ensembl.org/pub/current_{filetype}/{species}"
    
    url_strings = [base_url] + list(args)

    return "/".join(url_strings)


def read_file_line(filename):
    """Return the the first line of a file stripped of newlines.

    Args:
        filename (str): the path to the file
    """
    with open(filename) as f:
        return f.readline().strip()


def get_all_inputs():
    input_list = []

    input_list.extend(
        multiext(
            "resources/ensembl/genome",
            ".fa.gz",
            ".gtf.gz",
        )
    )

    # Files that apply to any organism
    input_list.append("resources/ensembl/star_genome")
    input_list.append("resources/ensembl/genome_info.json")
    input_list.append("resources/ensembl/gffutils.db")
    input_list.append("resources/ensembl/transcript_annotation.bed")

    # Files that only apply to human
    if config["ref"]["species"] == "homo_sapiens":
        input_list.append("resources/ensembl/star_genome_mane")
        # input_list.append("resources/ensembl/clinically_associated_variants.vcf.gz")
        input_list.append("resources/ensembl/cds_regions.tsv")
        input_list.append("resources/clinvar/clinvar.vcf.gz")
        input_list.append("resources/clinvar/clinvar.filtered.vcf.gz")
        input_list.append("resources/ensembl/MANE.GRCh38.v1.0.gtf.gz")

    return input_list