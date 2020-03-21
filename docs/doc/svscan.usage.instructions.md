# svscan tutorial

## How to run svscan with svscan

### 1. Prepare input bam

   <p>If your library is DNA, align your DNA library to hg19 with bwa mem, the hg19 fasta file can be any one which indexed by bwa mem.</p>
  
   <p>If your library is RNA, align your RNA library to refined hg19 transcriptome(```/share/work1/wulj/database/svdb/rna/ref/ref.fa```), other transcriptome is not supported.</p>

   <p>Your bwa must support SA tag operations. The alignment result BAM file must be sorted by coordinates, you can also mark duplications or remove duplications in the BAM file, svscan will not use any BAM record marked as duplication in furthur computation.</p>

### 2. Execute svscan and generate result file

### which svscan to use

the svscan ```/share/work1/wulj/programs/bin/svscan``` is for my personal development usage, which might not be stable for production usage, if you want to experience the new features or latest improvement of svscan, try this one.

the svscan ```/share/work1/svpapp/gitlab/svpipe/bin/svscan``` is used for production, which should be stable and will not crash on execution, if you want to do some production sample analysis, try this one.

### general example
execute ```/share/work1/svpapp/gitlab/svpipe/bin/svscan -h``` for help and do whatever you can, each options or flag has detailed explanations which should be easily understood.

### dna exmaple(with arguments used in production for gene set 86)
   
   ``` 
   /share/work1/svpapp/gitlab/svpipe/bin/svscan
   --whiteminsrs 1
   --usualminsrs 3
   --whitemindps 2
   --usualmindps 3
   --whiteminsrr 3
   --usualminsrr 5
   --whitemindpr 3
   --usualmindpr 5
   --whiteminttr 2
   --usualminttr 3
   --min_seed_sr 1
   --min_seed_dp 2
   -b /share/A30050T/OncDir/vannul/work/svpipe1/guoyushuai712/Z19W05624-B1TA/output/vori/06.dedup/pre_Z19W05624-B1TA_C58.mkdup.sort.bam
   -g /share/work1/wulj/database/hg19/hg19.fa
   --whitelist /share/work1/wulj/database/svdb/dna/wlist/fwlist.86.tsv
   --samegenel /share/work1/wulj/database/svdb/dna/slist/fslist.86.tsv
   --blacklist /share/work1/wulj/database/svdb/dna/blist/fblack.all.tsv
   --bgbcf /share/work1/wulj/database/svdb/dna/bgbcf/bgbcf.bcf
   --fusionrpt /share/A30050T/OncDir/vannul/work/svpipe1/guoyushuai712/Z19W05624-B1TA/output/vori/08.sv/pre_Z19W05624-B1TA_C58.fs.tsv
   --genecrdlist /share/work1/wulj/database/svdb/extra/gene.coord.tsv
   -c /share/work1/wulj/database/svdb/dna/regs/fbcreg.86.bed
   --extraanno /share/work1/wulj/database/svdb/dna/extra/fext.86.bed
   -a /share/work1/wulj/database/svdb/dna/anndb/anno.sort.gz
   -o /share/A30050T/OncDir/vannul/work/svpipe1/guoyushuai712/Z19W05624-B1TA/output/vori/08.sv/pre_Z19W05624-B1TA_C58.sv.bcf
   -t /share/A30050T/OncDir/vannul/work/svpipe1/guoyushuai712/Z19W05624-B1TA/output/vori/08.sv/pre_Z19W05624-B1TA_C58.sv.tsv
   -v /share/A30050T/OncDir/vannul/work/svpipe1/guoyushuai712/Z19W05624-B1TA/output/vori/08.sv/pre_Z19W05624-B1TA_C58.sv.bam
   -e /share/A30050T/OncDir/vannul/work/svpipe1/guoyushuai712/Z19W05624-B1TA/output/vori/08.sv/pre_Z19W05624-B1TA_C58.fr.xlsx
   -f /share/A30050T/OncDir/vannul/work/svpipe1/guoyushuai712/Z19W05624-B1TA/output/vori/08.sv/pre_Z19W05624-B1TA_C58.fr.tsv
   -n 8
   ```
  
###  rna example(with arguments used in production for gene set 105)
   
   ```
   /share/work1/svpapp/gitlab/svpipe/bin/svscan
   --whiteminsrs 1
   --usualminsrs 3
   --whitemindps 2
   --usualmindps 3
   --whiteminsrr 3
   --usualminsrr 5
   --whitemindpr 3
   --usualmindpr 5
   --whiteminttr 2
   --usualminttr 3
   --min_seed_sr 1
   --min_seed_dp 2
   --whitemindep 2
   --usualmindep 2
   --rna
   --idpdropmask 38003584
   --ndbdropmask 40101760
   -b /share/A30050T/OncDir/svpapp/wangx9211/WX0306/output/vori/06.dedup/PYZB305R4000_B28.mkdup.sort.bam
   -g /share/work1/svpapp/database/svdb/rna/ref/ref.fa
   --whitelist /share/work1/svpapp/database/svdb/rna/wlist/fwlist.105.tsv
   --samegenel /share/work1/svpapp/database/svdb/rna/slist/fslist.105.tsv
   --blacklist /share/work1/svpapp/database/svdb/rna/blist/fblack.all.tsv
   --fusionrpt /share/A30050T/OncDir/svpapp/wangx9211/WX0306/output/vori/08.sv/PYZB305R4000_B28.fs.tsv
   --genecrdlist /share/work1/svpapp/database/svdb/extra/gene.coord.tsv
   -c /share/work1/svpapp/database/svdb/rna/regs/fbcreg.105.bed
   --extraanno /share/work1/svpapp/database/svdb/rna/extra/fext.105.bed
   -a /share/work1/svpapp/database/svdb/rna/anndb/anno.sort.gz
   -o /share/A30050T/OncDir/svpapp/wangx9211/WX0306/output/vori/08.sv/PYZB305R4000_B28.sv.bcf
   -t /share/A30050T/OncDir/svpapp/wangx9211/WX0306/output/vori/08.sv/PYZB305R4000_B28.sv.tsv
   -v /share/A30050T/OncDir/svpapp/wangx9211/WX0306/output/vori/08.sv/PYZB305R4000_B28.sv.bam
   -e /share/A30050T/OncDir/svpapp/wangx9211/WX0306/output/vori/08.sv/PYZB305R4000_B28.fr.xlsx
   -f /share/A30050T/OncDir/svpapp/wangx9211/WX0306/output/vori/08.sv/PYZB305R4000_B28.fr.tsv
   -n 8
   ```
### database files
You may need to do analysis with some other gene set of DNA or RNA library, all the DNA and RNA database files available are in ```/share/work1/wulj/database/svdb```, they are distributed in a standardized directory tree and follow regulary nomenclature rules. In fact, you just replace the geneset names with the other ones, such as 105 to 31 in ```RNA example``` or 457 to 655 in ```DNA example``` will generate the commands needed to do 31 gene set RNA svscan and 655 gene set svscan anslysis.

All change history of database file can be tracked in its gitlab repo at [svdb](http://10.100.35.200:10080/wulj3253/svdb)

### generate svscan sh file automatically
Yeah, you can generate svscan bash script automatically. First put the bam paths in a plain text file with each path on one line, make sure no duplicated bam names appears in the same bamlist file.
#### automatically generate batch svscan scripts for a DNA bamlist
```/share/work1/wulj/programs/bin/prepDSV <in.list> <gset> <outdir> <optargs>```
#### automatically generate batch svscan scripts for a RNA bamlist
```/share/work1/wulj/programs/bin/prepRSV <in.list> <gset> <outdir> <optargs>```

|positional arguments| explanations
|--------------------|--------------
|in.list|you input bamlist
|gset|gene set of your bams in the input bamlist
|outdir|where to output all the scripts of you svscan
|optargs|optional arguments(does not needed)

after execution, you will get all the scripts in the outdir, and you can execute theme or qsub to computation nodes to run. all the arguments used are the same as svscan used in production. make sure the right version svscan is in your environment variable, if not, export the proper one before execution each script.


## How to run svscan with svpapp

Get svpapp of for your system from shared files in QQ group```和瑞研发生信-北京-HYS-```, or download them from server path ```/share/work1/wulj/database/release``` and install it.

Open svpapp->Help->Manual for help, you are welcome to run svscan with svapapp, which will save you a lot of valuable lifetime.
