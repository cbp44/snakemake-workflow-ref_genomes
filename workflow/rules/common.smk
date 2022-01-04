def get_ensembl_base_url(filetype, postfix=""):
    ensembl_release=config["ref"]["ensembl_release"]
    species=config["ref"]["species"]
    
    base_url = "https://ftp.ensembl.org/pub/release-{0}/{1}/{2}{3}".format(
        ensembl_release, filetype, species, "" if not postfix else "/"
    )

    base_url += postfix

    return base_url
