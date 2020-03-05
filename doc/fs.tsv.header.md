|Column| explanations
|------|-------------
SampleID|样品名
FusionGene|融合基因hgene->tgene(如果hgene或者tgene为-，则代表非基因区）
FusionReads|支持融合的分子（一个模版测到的两条读段只算一个分子）
TotalReads|支持融合分子数和参考型分子数总和
FusionRate|融合比率（FusionReads/TotalReads)
inDB|融合形式(hgene->tgene)是否在ONCOKB或者COSMIC中（Y表示在，N表示不在）
status|"异常状态marker<br>Y表示无异常<br>M表示tgene->hgene在ONCOKB或者COSMICD<br>D表示hgene和tgene不是正常的转录本水平hgene的5'部分正链和tgene3'部分正链相连<br>C表示hgene或者tgene中的某一个探针范围内基因不在其常见融合方向上<br>S表示hgene和tgene是相同的基因<br>G表示hgene 和tgene中有一个为非基因区"
fpDist|"基因间距离<br>0 表示hgene和tgene为相同基因，或者不同染色体上的基因或同一染色体上但是不相邻（中间有其他基因<br>负数表示hgene和tgene在基因组上有overlap<br>正数表示俩基因在基因组上相邻的距离，俩基因间没有其他基因"
Gene1|hgene名字
Chr1|hgene所在染色体名字
JunctionPosition1|hgene所在染色体上断点位置
Strand1|hgene所在染色体上的链向
Transcript1|hgene的转录本信息，格式为：转录本号，基因组上链向，单元属性，单元编号，外显子号。（单元是指内含子或者外显子）
Gene2|tgene名字
Chr2|tgene所在染色体名字
JunctionPosition2|tgene所在染色体上断点位置
Strand2|tgene所在染色体上的链向
Transcript2|tgene的转录本信息，格式为：转录本号，基因组上链向，单元属性，单元编号，外显子号。（单元是指内含子或者外显子）
FusionSequence|融合保守序列
fseqBp|融合保守序列断点位置
svType|"融合对应的结构变异事件类型<br>BND：不同染色体之间的结构变异<br>INV：倒位<br>INS：插入<br>|DEL：缺失<br>DUP：重复"
svSize|断点之间的距离，不同染色体之间的断点之间的距离为-1
srCount|断裂比对种子数（reads条数）
dpCount|不一致不对读段对种子数（分子数）
srRescued|所有支持融合的断裂比对reads数
dpRescued|所有支持融合的不一致比对读段对数（分子数）
srRefCount|支持参考型的断裂读段数（reads数）
dpRefCount|支持参考型的不一致比对读段对数（分子数）
insBp|断点附近插入序列长度
insSeq|断点附近插入序列
svID|SV ID
svtInt|结构变异事件分类0-8的整数:
fsMask|融合mask
fsHits|保守序列重比对返回值
ts1Name|hgene转录本名字（RNA融合特有）
ts1Pos|hgene转录本断点位置（RNA融合特有）
ts2Name|tgene转录本名字（RNA融合特有）
ts2Pos|tgene转录本断点位置（RNA融合特有）
fsCigar|断点两侧简单CIGAR值（RNA融合特有）
