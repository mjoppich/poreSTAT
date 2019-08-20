import HTSeq
import argparse
import datetime


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--input', '-i', type=argparse.FileType('rb'), required=True)
    parser.add_argument('--output', '-o', type=argparse.FileType('w'), required=True)
    args = parser.parse_args()


    almnt_file = HTSeq.BAM_Reader( args.input )

    counts =  HTSeq.GenomicArray( "auto", stranded=False,  typecode='i' )
    fcounts = HTSeq.GenomicArray( "auto", stranded=False,  typecode='i' )


    curChrom = None
    for almnt in almnt_file:
        if not almnt.aligned or almnt.not_primary_alignment or almnt.supplementary:
            continue

        if curChrom != almnt.iv.chrom:

            dt = datetime.datetime.now()
            
            print(dt.isoformat(), args.input.name, "Switching Chromosome", curChrom, almnt.iv.chrom)
            curChrom = almnt.iv.chrom

        for cigop in almnt.cigar:
            if not cigop.type == "M":
                continue

            counts[cigop.ref_iv] += 1

    print("Filtering counts")
    for iv, cnt in counts.steps():

        if cnt >= 10:
            fcounts[iv] = cnt

    print("Writing File")
    fcounts.write_bedgraph_file(args.output)

    #python3 createExpressedRegions.py --input /mnt/f/schober_data/Mlet7/Mlet7A1/star_sorted.bam --output /mnt/f/schober_data/Mlet7/Mlet7A1/counts.wig