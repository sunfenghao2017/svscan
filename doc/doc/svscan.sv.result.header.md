|Column|解释|Explanation
|------|--------------|--------------
svType|融合对应的结构变异事件类型<br>BND：不同染色体之间的结构变异<br>INV：倒位<br>INS：插入<br>DEL：缺失<br>DUP:重复<br>|structural variant event type of this fusion<br>BND: translocation across chrosome/contig<br>INV: inversion<br>INS: insertion<br>DEL: deletion<br>DUP: duplication
svSize|断点之间的距离，不同染色体之间的断点之间的距离为-1|breakpoint distance(-1 if breakpoints are on different chromosome)
bpMark|断点标记<br>L:倒位左侧断点<br>R:倒位右侧断点<br>5'->5':不同Ref之间5‘端和5’端连结型断点<br>5'->3':不同Ref之间5‘端和3’端连结型断点<br>3'->5':不同Ref之间3‘端和5’端连结型断点<br>3'->3':不同Ref之间3‘端和3’端连结型断点|breakpoint marker<br>L: leftmost breakpoint of inversion<br>R: rightmost breakpoint of inversion<br>5'->5':5' 5' part of two part of contig connected<br>5'->3': 5'part and 3'part of contig connected<br>3'->5': 3'part and 5'part of contig connected<br>3'->3': 3' part of two part of contig connected
bp1Chr|第一个断点所在的染色体|chromosome breakpoint1 is on
bp1Pos|第一个断点位置|breakpoint position of breakpoint1
bp2Chr|第二个断点所在的染色体|chromosome breakpoint1 is on
bp2Pos|第二个断点位置|breakpoint position of breakpoint2
srCount|断裂比对种子数(reads条数)|split read seed count(in read)
dpCount|不一致不对读段对种子数(分子数)|discordant pair seed count(in molecule)
srRescued|所有支持融合的断裂比对reads数|all split reads supporting this fusion(in read)
dpRescued|所有支持融合的不一致比对读段对数(分子数)|all discordant pairs supporting this fusion(in molecule)
srRefCount|支持参考型的断裂读段数(reads数)|all split reads supportint reference type(in read)
dpRefCount|支持参考型的不一致比对读段对数(分子数)|all discordant pairs supporting this fusion(in molecule)
molRescued|支持融合(结构变异)的分子(一个模版测到的两条读段只算一个分子)|all molecules support this sv(fusion) event
AF|融合(结构变异)比率(molRescued/(max(srRefCount,  dpRefCount)+ molRescued)|rate = （molRescued/(max(srRefCount, dpRefCount)+ molRescued)
srSRescued|拯救回来的断裂比对种子数(reads条数)|rescued split read seed count(in read)
tsrSResMaln|拯救回来的断裂比对种子中探针设计基因的partner基因发生多位置比对的数量(reads条数)|number of  rescued split read seeds with its split patner(which is not on the probe capture region mapping) mapping on at least 2 different position on reference
srSResMalnRate|tsrSResMaln/srSRescued|tsrSResMaln/srSRescued
insBp|断点附近插入序列长度|length of sequence inserted after breakpoint
insSeq|断点附近插入序列|insertion sequence after breakpoint
svSeq|融合(结构变异)保守序列|consensus sequence from MSA of seed split reads supporting this fusion
seqBp|融合(结构变异)保守序列断点位置|breakpoint position on the consensus sequence when aligned to reference
ID|SV ID|SV ID(for internal and across table reference usage)
svtInt|结构变异事件分类0-8的整数|integer representation of the catenation of two translocated segment of two DNA molecules<br>0: intra contig 5'->5'<br>1: intra contig 3'->3'<br>2: intra contig 5'->3'<br>3: intra contig 3'->5'<br>4: no translocation<br>5: inter contig 5'->5'<br>6: inter contig 3'->3'<br>7: inter contig 5'->3'<br>8: inter contig 3'->5'bp1Gene|第一个断点位置上的基因信息
bp1Gene|第一个断点位置上的基因信息|gene information of breakpoint1
bp2Gene|第二个断点位置上的基因信息|gene information of breakpoint2
fuseGene|该结构变异事件导致的融合基因信息|fusion gene information this event supports
fsMask|融合mask|fusion property bit mask value
fsHits|保守序列重比对返回值|consensus sequence of this fusion event realignment return value<br>0: nice perfect match or just one primary sc<br>-1: whole seq match continuous on reference<br>-2: two primary compatible softclip match<br>-3: breakpoint conflicts from original breakpoint<br>>=4: number of compatible split alignment on reference
ts1Name|hgene转录本名字(RNA融合特有)|transcript name of hgene(RNA mode only)
ts1Pos|hgene转录本断点位置(RNA融合特有)|breakpoint position of hgene transcript(RNA mode only)
ts2Name|tgene转录本名字(RNA融合特有)|transcript name of tgene(RNA mode only)
ts2Pos|tgene转录本断点位置(RNA融合特有)|breakpoint position of tgene transcript(RNA mode only)
fsCigar|断点两侧简单CIGAR值(RNA融合特有)|consensus sequence alignment CIGAR around breakpoint(RNA mode only)
