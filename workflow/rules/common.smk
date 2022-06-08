def get_ensembl_url(filetype, *args):
    """Get a full Ensembl download URL for the configured species given a
    specific filetype (like fasta, gtf, etc.) and optional postfix strings.
    """
    species = config["ref"]["species"]
    base_url = f"https://ftp.ensembl.org/pub/current_{filetype}/{species}"
    
    url_strings = [base_url] + list(args)

    return "/".join(url_strings)


def read_file_line(filename):
    """Return the the first line of a file stripped of newlines.

    Args:
        filename (str): the path to the file
    """
    with open(filename) as f:
        return f.readline().strip()