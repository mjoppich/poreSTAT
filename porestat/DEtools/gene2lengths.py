#!/usr/bin/env python3
"""
gene2lengths.py
Markus Joppich
August 2019

This work uses two python-scripts from Kamil Slowikowski https://slowkow.com/notes/gencode-coding-lengths/ which are licenced under CC4 BY-SA .
In the same fashion, this adaptation is licenced under CC4 BY-SA: https://creativecommons.org/licenses/by-sa/4.0/

I have made this code python3 compatible and added argparse .

For my configuration I did not need the gene2ensembl conversion, hence I commented it out. I added a gene_biotype output in the final output as I thought that'd be helpful in sorting out rRNA genes.


coding_lengths.py
Kamil Slowikowski
February 7, 2014
Count the number of coding base pairs in each Gencode gene.
Gencode coordinates, including all exons with Ensembl identifiers.
(Gencode release 17 corresponds to hg19)
    ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
    ftp://ftp.sanger.ac.uk/pub/gencode/release_17/gencode.v17.annotation.gtf.gz
    chr1  HAVANA  gene        11869  14412  .  +  .  gene_id "ENSG00000223972.4";
    chr1  HAVANA  transcript  11869  14409  .  +  .  gene_id "ENSG00000223972.4";
    chr1  HAVANA  exon        11869  12227  .  +  .  gene_id "ENSG00000223972.4";
    chr1  HAVANA  exon        12613  12721  .  +  .  gene_id "ENSG00000223972.4";
    chr1  HAVANA  exon        13221  14409  .  +  .  gene_id "ENSG00000223972.4";
NCBI mapping from Entrez GeneID to Ensembl identifiers.
    ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
    9606  1  ENSG00000121410  NM_130786.3     ENST00000263100  NP_570602.2     ENSP00000263100
    9606  2  ENSG00000175899  NM_000014.4     ENST00000318602  NP_000005.2     ENSP00000323929
    9606  3  ENSG00000256069  NR_040112.1     ENST00000543404  -               -
    9606  9  ENSG00000171428  NM_000662.5     ENST00000307719  NP_000653.3     ENSP00000307218
    9606  9  ENSG00000171428  XM_005273679.1  ENST00000517492  XP_005273736.1  ENSP00000429407
Output:
    Ensembl_gene_identifier  GeneID  length
    ENSG00000000005          64102   1339
    ENSG00000000419          8813    1185
    ENSG00000000457          57147   3755
    ENSG00000000938          2268    3167
USAGE:
    gene2lengths.py -g FILE -n FILE [-o FILE]
OPTIONS:
    -h          Show this help message.
    -g FILE     Gencode annotation.gtf.gz file.
    -n FILE     NCBI gene2ensembl.gz file.
    -o FILE     Output file (gzipped).
"""


"""
GTF.py
Kamil Slowikowski
December 24, 2013

Read GFF/GTF files. Works with gzip compressed files and pandas.

    http://useast.ensembl.org/info/website/upload/gff.html

LICENSE

This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
"""


from collections import defaultdict
import gzip
import pandas as pd
import re


GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')

class GTF:

    @classmethod
    def dataframe(cls, filename):
        """Open an optionally gzipped GTF file and return a pandas.DataFrame.
        """
        # Each column is a list stored as a value in this dict.
        result = defaultdict(list)

        for i, line in enumerate(cls.lines(filename)):
            for key in line.keys():
                # This key has not been seen yet, so set it to None for all
                # previous lines.
                if key not in result:
                    result[key] = [None] * i

            # Ensure this row has some value for each column.
            for key in result.keys():
                result[key].append(line.get(key, None))

        return pd.DataFrame(result)

    @classmethod
    def lines(cls, filename):
        """Open an optionally gzipped GTF file and generate a dict for each line.
        """
        fn_open = gzip.open if filename.endswith('.gz') else open

        with fn_open(filename) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                else:
                    yield cls.parse(line)

    @classmethod
    def parse(cls, line):
        """Parse a single GTF line and return a dict.
        """
        result = {}

        fields = line.rstrip().split('\t')

        for i, col in enumerate(GTF_HEADER):
            result[col] = cls._get_value(fields[i])

        # INFO field consists of "key1=value;key2=value;...".
        infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

        for i, info in enumerate(infos, 1):
            # It should be key="value".
            try:
                key, _, value = re.split(R_KEYVALUE, info, 1)
            # But sometimes it is just "value".
            except ValueError:
                key = 'INFO{}'.format(i)
                value = info
            # Ignore the field if there is no value.
            if value:
                result[key] = cls._get_value(value)

        return result

    @classmethod
    def _get_value(cls, value):
        if not value:
            return None

        # Strip double and single quotes.
        value = value.strip('"\'')

        # Return a list if the value has a comma.
        if ',' in value:
            value = re.split(R_COMMA, value)
        # These values are equivalent to None.
        elif value in ['', '.', 'NA']:
            return None

        return value



import argparse
import pandas as pd
import gzip
import time
import sys
from contextlib import contextmanager


def main(args):
    # Input files.
    GENCODE      = args.gtf.name
    NCBI_ENSEMBL = None

    # Output file prefix.
    GENE_LENGTHS = args.output.name#"ncbi_ensembl_coding_lengths.tsv"

    with log("Reading the Gencode annotation file: {}".format(GENCODE)):
        gc = GTF.dataframe(GENCODE)

    # Select just exons of protein coding genes, and columns that we want to use.
    idx = (gc.feature == 'exon')# & (gc.transcript_biotype == 'protein_coding')

    exonData = ['seqname','start','end','gene_id']

    if "gene_biotype" in gc:
        exonData.append("gene_biotype")

    exon = gc.ix[idx, exonData] #gene_name

    # Convert columns to proper types.
    exon.start = exon.start.astype(int)
    exon.end = exon.end.astype(int)

    # Sort in place.
    exon.sort_values(by=['seqname','start','end'], inplace=True)

    # Group the rows by the Ensembl gene identifier (with version numbers.)
    groups = exon.groupby('gene_id')

    with log("Calculating coding region (exonic) length for each gene..."):
        lengths = groups.apply(count_bp)
        bioType = groups.apply(get_biotype)

    if NCBI_ENSEMBL != None:
        with log("Reading NCBI mapping of Entrez GeneID to Ensembl gene identifier: {}".format(NCBI_ENSEMBL)):
            g2e = pd.read_table(NCBI_ENSEMBL,
                                compression="gzip",
                                header=None,
                                names=['tax_id', 'GeneID',
                                    'Ensembl_gene_identifier',
                                    'RNA_nucleotide_accession.version',
                                    'Ensembl_rna_identifier',
                                    'protein_accession.version',
                                    'Ensembl_protein_identifier'])

    # Create a new DataFrame with gene lengths and EnsemblID.
    ensembl_no_version = lengths.index.map(lambda x: x.split(".")[0])
    ldf = pd.DataFrame({'length': lengths,
                        'gene_identifier': ensembl_no_version,
                        'biotype': bioType},
                        index=lengths.index)

    # Merge so we have EntrezGeneID with length.
    #m1 = pd.merge(ldf, g2e, on='Ensembl_gene_identifier')
    #m1 = m1[['Ensembl_gene_identifier', 'GeneID', 'length']].drop_duplicates()
    m1 = ldf.drop_duplicates()

    with log("Writing output file: {}".format(GENE_LENGTHS)):
        with open(GENE_LENGTHS, "w") as out:
            m1.to_csv(out, sep="\t", index=False)


def get_biotype(df):

    if not "gene_biotype" in df:
        return ""

    return "|".join(set(df["gene_biotype"]))


def count_bp(df):
    """Given a DataFrame with the exon coordinates from Gencode for a single
    gene, return the total number of coding bases in that gene.
    Example:
        >>> import numpy as np
        >>> n = 3
        >>> r = lambda x: np.random.sample(x) * 10
        >>> d = pd.DataFrame([np.sort([a,b]) for a,b in zip(r(n), r(n))], columns=['start','end']).astype(int)
        >>> d
           start  end
        0      6    9
        1      3    4
        2      4    9
        >>> count_bp(d)
        7
    Here is a visual representation of the 3 exons and the way they are added:
          123456789  Length
        0      ----       4
        1   --            2
        2    ------       6
            =======       7
    """
    start = df.start.min()
    end = df.end.max()
    bp = [False] * (end - start + 1)
    for i in range(df.shape[0]):
        s = df.iloc[i]['start'] - start
        e = df.iloc[i]['end'] - start + 1
        bp[s:e] = [True] * (e - s)
    return sum(bp)


@contextmanager
def log(message):
    """Log a timestamp, a message, and the elapsed time to stderr."""
    start = time.time()
    sys.stderr.write("{} # {}\n".format(time.asctime(), message))
    yield
    elapsed = int(time.time() - start + 0.5)
    sys.stderr.write("{} # done in {} s\n".format(time.asctime(), elapsed))
    sys.stderr.flush()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some gtf.')
    parser.add_argument("-g", "--gtf", type=argparse.FileType("r"), required=True)
    parser.add_argument("-o", "--output", type=argparse.FileType("w"), required=True)
    args = parser.parse_args()


    main(args)