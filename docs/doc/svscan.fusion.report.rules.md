# svscan fusion report rules

```appliable to svscan 0.1.11-r5 and later version```

```The following rules will take the fusion hgene->tgene as an example.```

1. if hgene or tgene is in the [gene blacklist](http://10.100.35.200:10080/wulj3253/svdb/tree/master/dna/blist/fblack.all.tsv), it will not be reported
2. if hgene->tgene is in the [fusion blacklist](http://10.100.35.200:10080/wulj3253/svdb/tree/master/dna/blist/fblack.all.tsv), it will not be reported
3. if hgene->tgene is in the [background pool](http://10.100.35.200:10080/wulj3253/svdb/tree/master/dna/bgbcf/bgbcf.bcf), it will not be reported
4.	if hgene->tgene is identified from split alignment reads and hgene->tgene is not in [database](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv) and the consensus sequence around breakpoint is some simple repeat sequence is it will not be reported.<br>any sequende met any one of the following two conditions will be defined as simple repeat sequence<br>1) count two base count of ```ATCG``` is greater than 0.85 * total length of the sequence<br>2) the sequence contains at least 16 consecutive ```GT```, ```TG```, ```AC``` or ```CA```
5.	if hgene tgene are the same gene, and the structural variant type of  hgene->tgene is not in the predefined intra-gene report range([DNA](http://10.100.35.200:10080/wulj3253/svdb/tree/master/dna/slist), [RNA](http://10.100.35.200:10080/wulj3253/svdb/tree/master/rna/slist),) it will not be reported
6. if structural variant size is too small, it will not be reported<br> structural variant size is too small if one of the following condidions satisfied:<br>1) structural variant is in the same intron of the same gene<br>2) structural variant is in the same exon of the same gene but the distance between two breakpoint is less than 0.8 of the exon length<br>3) if one partner in hgene->tgene is on same chromosome and the distance between two breakpoint is less than 2k
7. if the fusion is not in [database](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv) and its af is less than 0.01(tissue) or less than 0.005(plasma and hydrothorax), it will not be reported, however i'm considering report some low frequency fusion with ```really good supporting reads``` which mets the following conditions:<br>1) The fusion gene transcript is functional<br>2) Split reads count is 2times of the normal seed requirement<br>3) Rescued split reads count is 2times of the normal seed requirement<br>4) Supporting molecules is 1.5times of the normal support molecules
8. threshold of breakpoint depth: fusion in [database](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv) must be greater than 30x, other fusions must be greater than 300x(does not applies to RNA fusion)
9. threshold of supporting molecules: fusion in [database](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv)  must be greater than 2, other fusions must be greater than 3
10. threshold of seed count: fusion in [database](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv) must has at least 3 split read support or 2 discordant pairr support, other fusions must has at least 3 split read support or 3 discordant pair support
11. threshold of total read count: fusion in [database](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv)  must has at least 3 split read support or 2 discordant pairs, other fusions must has at least 5 split read support or 5 discordant pair support
12. if hgene->tgene is not in the predefined report range, it will not be reported
13. at least one of hgene and tgene is a gene and in the panel probe range
14. if hgene->tgene is not in [database](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv) and the consensus sequence does not map well, it will not be reported.<br>if any one of the following 4 conditions satisfied, the consensus sequene mapping is not good<br>1) the consensus sequence is mapped on the reference consecutively without any split<br>2) there are more than 2 pairs of compatible split alignment<br>3) there are more than 2 pairs of primary split alignment<br>4) the only one primary split alignment breakpoint positions conflict with the original positions
16. if there are multiple pattern of hgene->tgene fusion in on sample, the most reasonable one will be reported
17. NCRNA gene participated fusion will not be reported in RNA mode
18. if the target gene partner part of rescued split reads are all in repeat region, it will not be reported
