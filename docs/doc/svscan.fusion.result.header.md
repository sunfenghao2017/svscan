|Column|解释|Explanation
|------|--------------|-------------
FusionGene|融合基因hgene->tgene(如果hgene或者tgene为-，则代表非基因区)|fusion gene in the format hgene->tgene(if one of the partner is '-', then it stands for a non-gene region)
FusionMols|支持融合的分子(一个模版测到的两条读段只算一个分子)|number of molecules supporting fusion(two paired read from the same template should only be counted as one molecule
TotalMols|支持融合分子数和参考型分子数总和|number of molecules supporting fusion and reference type
FusionRate|融合比率（FusionMols/TotalMols)|fusion rate, equals FusionReads/TotalReads
inDB|融合形式hgene->tgene或tgene->hgene是否在[数据库](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv)（Y表示在，N表示不在)|whether the hgene->tgene(tgene->hgene) fusion pattern exists in [database](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv)(Y for positive, N for negative)
fsPattern|如果是断裂比对read支持的融合，第一个符号表示支持该融合的read比对到第一个基因所在基因组上的时候，比对到的基因组链向，第二个符号表示该read比对到第二个基因所在基因组上的时候，比对到的基因组的链向<br>如果是不一致比对读段对支持的融合，第一个符号表示支持改融合的一条read比对到第一个基因所在基因组上的时候，比对到的基因组链向，第二个符号表示另一条read比对到第二个基因所在基因组上的时候，比对到基因组的链向|sign to stand for the read alignenment strand on genome of each fusion partner, first sign stands for the alignment strand on genome of hgene(SR or DP), second sign stands for the alignment strand on genome of tgene(SR or DP)
status|异常状态marker<br>Y表示无异常<br>M表示tgene->hgene在[database](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv)<br>D表示hgene和tgene不是正常的转录本水平hgene的5'部分正链和tgene3'部分正链相连<br>C表示hgene或者tgene中的某一个探针范围内基因不在其常见融合方向上<br>S表示hgene和tgene是相同的基因<br>G表示hgene 和tgene中有一个为非基因区|abnormal status marker<br>Y means normal<br>M means tgene->hgene exists in [database](http://10.100.35.200:10080/wulj3253/svdb/blob/master/extra/fusedb.tsv)<br>D means hgene->tgene is not a standard fusion transcript(a standard transcript transcript of hgene->tgene should consists of 5' sense strand of hgene and 3' sense strand of tgene)<br>C means either tgene or hgene appears not in the common direction of reported fusion events that gene participated<br>S means hgene and tgene is the same gene<br>G means at least one of hgene and tgene is the non-gene region
fpDist|基因间距离<br>0 表示hgene和tgene为相同基因，或者不同染色体上的基因或同一染色体上但是不相邻（中间有其他基因<br>负数表示hgene和tgene在基因组上有overlap<br>正数表示俩基因在基因组上相邻的距离，俩基因间没有其他基因|distance of hgene and tgene<br>0 means hgene and tgene is the same gene or hgene and tgene are from different chromosome or they are not near each other even through they are on the same chromosome<br>negative value means hgene and tgene overlaps on the same chromosome<br>positive value means that tgene and hgene are near each other without any gene between them
Gene1|hgene名字|RefGene name of hgene
Chr1|hgene所在染色体名字|chrosome name of hgene
JunctionPosition1|hgene所在染色体上断点位置(```1 based```)|breakpoint position of hgene(```1 based```)
Strand1|hgene所在染色体上的链向|strand of chromosome hgene is on
Transcript1|hgene的转录本信息，格式为：转录本号，基因组上链向，单元属性，单元编号，外显子号.(单元是指内含子或者外显子)|transcript information of hgene in the format ```transcript number,strand,unit name,unit number,exon number```(uint is just exon or intron)
Gene2|tgene名字|RefGene name of tgene
Chr2|tgene所在染色体名字|chrosome name of tgene
JunctionPosition2|tgene所在染色体上断点位置(```1 based```)|breakpoint position of tgene(```1 based```)
Strand2|tgene所在染色体上的链向|strand of chromosome tgene is on
Transcript2|tgene的转录本信息，格式为：转录本号，基因组上链向，单元属性，单元编号，外显子号.(单元是指内含子或者外显子)|transcript information of tgene in the format ```transcript number,strand,unit name,unit number,exon number```(uint is just exon or intron)
FusionSequence|融合保守序列|consensus sequence from MSA of seed split reads supporting this fusion
seqBp|融合(结构变异)保守序列断点位置(```1 based```)，该位置后面紧跟着的序列（为插入序列，如果有)和另一个partner基因的序列|breakpoint position on the consensus sequence when aligned to reference(```1 based```), after the position, the sequence belongs to the other partner excluding the insertions sequence if there is
svType|融合对应的结构变异事件类型<br>BND：不同染色体之间的结构变异<br>INV：倒位<br>INS：插入<br>DEL：缺失<br>DUP:重复<br>|structural variant event type of this fusion<br>BND: translocation across chrosome/contig<br>INV: inversion<br>INS: insertion<br>DEL: deletion<br>DUP: duplication
svSize|断点之间的距离，不同染色体之间的断点之间的距离为-1|breakpoint distance(-1 if breakpoints are on different chromosome)
srCount|断裂比对种子数(reads条数)|split read seed count(in read)
dpCount|不一致比对读段对种子数(分子数)|discordant pair seed count(in molecule)
srRescued|所有支持融合的断裂比对reads数|all split reads supporting this fusion(in read)
dpRescued|所有支持融合的不一致比对读段对数(分子数)|all discordant pairs supporting this fusion(in molecule)
srRefCount|支持参考型的断裂读段数(reads数)|all split reads supportint reference type(in read)
dpRefCount|支持参考型的不一致比对读段对数(分子数)|all discordant pairs supporting this fusion(in molecule)
srSRescued|拯救回来的断裂比对种子数(reads条数)|rescued split read seed count(in read)
tsrSResMaln|拯救回来的断裂比对种子中探针设计基因的partner基因发生多位置比对的数量(reads条数)|number of  rescued split read seeds with its split patner(which is not on the probe capture region mapping) mapping on at least 2 different position on reference
srSResMalnRate|tsrSResMaln/srSRescued|tsrSResMaln/srSRescued
insBp|断点附近插入序列长度|length of sequence inserted after breakpoint
insSeq|断点附近插入序列|insertion sequence after breakpoint
svID|SV ID|SV ID(for internal and across table reference usage)
svtInt|结构变异事件分类0-8的整数|integer representation of the catenation of two translocated segment of two DNA molecules<br>0: intra contig 5'->5'<br>1: intra contig 3'->3'<br>2: intra contig 5'->3'<br>3: intra contig 3'->5'<br>4: no translocation<br>5: inter contig 5'->5'<br>6: inter contig 3'->3'<br>7: inter contig 5'->3'<br>8: inter contig 3'->5'
fsMask|融合mask|fusion property bit mask value
fsHits|保守序列重比对返回值|consensus sequence of this fusion event realignment return value<br>0: nice perfect match or just one primary sc<br>-1: whole seq match continuous on reference<br>-2: two primary compatible softclip match<br>-3: breakpoint conflicts from original breakpoint<br>>=4: number of compatible split alignment on reference
ts1Name|hgene转录本名字(RNA融合特有)|transcript name of hgene(RNA mode only)
ts1Pos|hgene转录本断点位置(RNA融合特有,```1 based```)|breakpoint position of hgene transcript(RNA mode only,```1 based```)
ts2Name|tgene转录本名字(RNA融合特有)|transcript name of tgene(RNA mode only)
ts2Pos|tgene转录本断点位置(RNA融合特有,```1 based```)|breakpoint position of tgene transcript(RNA mode only,```1 based```)
fsCigar|断点两侧简单CIGAR值(RNA融合特有)|consensus sequence alignment CIGAR around breakpoint(RNA mode only)
exon1|推测的hgene参与融合的外显子号，非基因区为-1（DNA融合特有）|estimated exon number of hgene participating the fusion event(DNA mode only)
exon2|推测的tgene参与融合的外显子号，非基因区为-1（DNA融合特有）|estimated exon number of tgene participating the fusion event(DNA mode only)