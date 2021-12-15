def get_ensembl_base_url(ensembl_release, filetype, species, postfix=""):
    base_url = "https://ftp.ensembl.org/pub/release-{0}/{1}/{2}{3}".format(
        ensembl_release,
        filetype,
        species,
        "" if not postfix else "/" + postfix
    )
    
    return base_url
