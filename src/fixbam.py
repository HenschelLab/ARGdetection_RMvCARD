import pysam
import sys

inbam, outbam = sys.argv[1:3]
"AFY01.sorted.length_100.bam"
sf=pysam.AlignmentFile(inbam, "rb")

## extracting the header, modifying it
for read in sf.fetch():
        hd = read.header.as_dict()
        for refgen in hd['SQ']:
            refgen['SN'] = refgen['SN'].replace('(',"_").replace(')',"_").replace("'", "_").replace(",", "_")        
        break

sf.close()
sf=pysam.AlignmentFile(inbam, "rb")
ff=pysam.AlignmentFile(outbam, "wb", header=hd)
for read in sf.fetch():
    a = pysam.AlignedSegment(ff.header)
    a.reference_name = read.reference_name.replace('(',"_").replace(')',"_").replace("'", "_").replace(",", "_")
    a.query_name = read.query_name
    a.query_sequence = read.query_sequence
    a.flag = read.flag
    a.reference_start = read.reference_start
    a.mapping_quality = read.mapping_quality
    a.cigar = read.cigar
    a.next_reference_id = read.next_reference_id
    a.next_reference_start = read.next_reference_start
    a.template_length = read.template_length
    a.query_qualities = read.query_qualities
    a.tags = read.tags
    ff.write(a)
ff.close()
sf.close()


#    if ')' in read.reference_name:
#        
#        read.header.from_dict(hd)

#    ff.write(read)

#

