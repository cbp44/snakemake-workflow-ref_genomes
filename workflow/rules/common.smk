def get_ensembl_base_url(filetype, postfix=""):
    # ensembl_release = config["ref"]["ensembl_release"]
    species = config["ref"]["species"]
    base_url = "https://ftp.ensembl.org/pub/current_{0}/{1}/{2}{3}".format(
        filetype, species, postfix, "" if not postfix else "/"
    )

    # base_url += postfix

    return base_url


def read_file_line(filename):
    """Return the the first line of a file stripped of newlines.

    Args:
        filename (str): the path to the file
    """
    with open(filename) as f:
        return f.readline().strip()