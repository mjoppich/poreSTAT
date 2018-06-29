import os

class FASTQ:

    def __init__(self, id, seq, qual):

        self.seq = seq.strip()
        self.id = id.strip()

        if id[0] == '@':
            self.id = id[1:]

        self.qual = qual.strip()

    def __str__(self):
        return '@' + self.id + os.linesep + self.seq + os.linesep + '+' + os.linesep + self.qual

    def __len__(self):
        return len(self.seq)

    @classmethod
    def parseFromStr(cls, fqrec):

        afq = fqrec.split("\n")

        return FASTQ(afq[0], afq[1], afq[3])

    def to_fasta(self):

        return FASTA(self.id, self.seq)

class FASTA:

    def __init__(self, id, seq):

        self.seq = seq.strip()
        self.id = id.strip()

        if id[0] == '>':
            self.id = id[1:]

    def __str__(self):

        return '>' + self.id + os.linesep + self.seq

    def __len__(self):

        return len(self.seq)