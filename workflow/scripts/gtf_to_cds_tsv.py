import gffutils

db = gffutils.FeatureDB(snakemake.input.db)

with open(snakemake.output.tsv, 'w') as outfileobj:
    outfileobj.write('{}\n'.format('\t'.join(['chr','start','end','rank','ensembl','refseq','gene_hgnc'])))
    for cds in db.features_of_type("CDS", order_by="start"):
        bed = [s.strip() for s in db.bed12(cds).split('\t')]
        if "gene_id" in cds.attributes and "gene_name" in cds.attributes:
            out_table = [bed[0],bed[1],bed[2],"-",cds.attributes["gene_id"][0],"-",cds.attributes["gene_name"][0]]
            outfileobj.write('{}\n'.format('\t'.join(out_table)))
