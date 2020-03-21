
column|explanations
------|---------------
svid|id of structural variant event this read supports
fsgene|fusion gene
rmap|primary alignment status of this read(```chr,pos,strand,cigar```)
mmap|primary alignment status of mate read(```chr,pos,strand,cigar```)
sa|supplementary alignment status of this read(```chr,pos,strand,cigar```)
rbp|breakpoint identified by the primary alignment of this read(```-1 if this read does not support the structural variant event as split read```)
sbp|breakpoint identified by the supplementary alignment of this read(```-1 if this read does not support the structural variant event as split read```)
lhit|number of mapping places of leading sequence on genome(if not evaluated, it is 0)
thit|number of mapping places of tailing sequence on genome(if not evaluated, it is 0)
lseq|leading sequence before breakpoint of this read
tseq|tailing sequence after breakpoint of this read
barcode|barcode of this read if BC tag exists in the bamrecord else 0
qname|read name of this read
read1|true if this is read1, false if this is read2
svrt|0 if this read supports the structural variant event in split read mode; 1 if this read supports the structural variant event in discordant pair mode
