from math import log2

seqlen = 0
nbases = 0

with open(snakemake.input[0]) as fin:
    line = fin.readline().strip()
    print("line: {}".format(line))
    seqlen = int(line)

print("sequence length: {}".format(seqlen))

nbases = round(min(14, (log2(seqlen) / 2) - 1))

print("nbases: {}".format(nbases))

if nbases > 0:
    with open(snakemake.output[0], "w") as fout:
        fout.write("{}".format(nbases))