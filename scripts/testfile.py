from time import strftime, gmtime

import pysam
from Bio import SeqIO

refSeqs = {}

with open("/mnt/c/igem/ecoli/ecoli.phage.fasta", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):

        refSeqs[record.id] = str(record.seq)

samfile = pysam.AlignmentFile("/mnt/c/igem/ecoli/ecoli.phage.sam", "r")

stats = {}

for x in samfile.header.references:
    print(x, samfile.header.get_reference_length(x), x in refSeqs)

    cov = [0] * samfile.header.get_reference_length(x)
    mm = [0] * samfile.header.get_reference_length(x)

    stats[x] = (cov, mm)


readIdx = 0
totalIdx = 0
for aln in samfile:

    totalIdx += 1

    if aln.is_unmapped:
        continue

    readIdx += 1

    if readIdx % 10000 == 0:
        print(strftime("%H:%M:%S", gmtime()), "Scanning read", readIdx, "Total Count", totalIdx)

    refSeq = refSeqs[aln.reference_name]
    allPairs = aln.get_aligned_pairs(matches_only=True)

    ms = stats[aln.reference_name]

    for readPos, seqPos in allPairs:
        ms[0][seqPos] += 1

        if aln.seq[readPos] != refSeq[seqPos]:
            ms[1][seqPos] += 1


with open("/mnt/c/igem/ecoli/ecoli.phage.coverage.bed", 'w') as fout:

    for x in stats:
        cov = stats[x][0]
        for i in range(0, len(cov)):
            fout.write("\t".join([str(x), str(i), str(i+1), str(cov[i])]) + "\n")

with open("/mnt/c/igem/ecoli/ecoli.phage.bed", 'w') as fout:

    for x in stats:
        mm = stats[x][1]
        for i in range(0, len(mm)):
            fout.write("\t".join([str(x), str(i), str(i+1), str(mm[i])]) + "\n")