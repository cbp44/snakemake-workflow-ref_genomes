import gffutils

db = gffutils.FeatureDB(snakemake.input.db)

with open(snakemake.output.bed, 'w') as outfileobj:
    for tx in db.features_of_type('transcript', order_by='start'):
        bed = [s.strip() for s in db.bed12(tx).split('\t')]
        bed[3] = tx.id
        outfileobj.write('{}\n'.format('\t'.join(bed)))
