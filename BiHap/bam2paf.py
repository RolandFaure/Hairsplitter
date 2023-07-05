import pysam as ps
import sys

def convert_bam_to_paf(bam_file, paf_file):
    infile = ps.AlignmentFile(bam_file, 'rb')
    outfile = open(paf_file, 'w')
    for read in infile:
        strand = "+"
        if read.is_reverse:
            strand = "-"
        blockLength = read.query_alignment_end - read.query_alignment_start

        if not read.is_unmapped :
            outfile.write(read.query_name+ "\t"\
                            + str(read.infer_read_length())+ "\t"\
                            + str(read.query_alignment_start)+ "\t"\
                            + str(read.query_alignment_end)+ "\t" \
                            + strand + "\t"\
                            + read.reference_name+ "\t"\
                            + str(infile.header.get_reference_length(read.reference_name))+ "\t"\
                            + str(read.reference_start)+ "\t"\
                            + str(read.reference_end)+ "\t"\
                            + str(read.get_cigar_stats()[0][0])+"\t"\
                            + str(blockLength)+  "\t"\
                            + str(read.mapping_quality)+ "\n")
        


if __name__ == "__main__":
    #the first argument is the bam file, the second is the paf file
    bamfile = sys.argv[1]
    paffile = sys.argv[2]

    convert_bam_to_paf(bamfile, paffile)
