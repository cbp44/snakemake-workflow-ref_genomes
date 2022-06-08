import gffutils

db = gffutils.create_db(snakemake.input[0],
                        dbfn=snakemake.output[0],
                        force=True,
                        keep_order=False,
                        merge_strategy="warning",
                        sort_attribute_values=False,
                        disable_infer_genes=True,
                        disable_infer_transcripts=True)