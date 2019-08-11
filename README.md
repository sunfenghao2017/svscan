program: sver  
version: 0.0.0  
updated: 18:05:08 Aug 11 2019  
Usage: sver [OPTIONS]  

 
|Options                           | Explanations
|----------------------------------|-----------------------------------------------------
|  -h,--help                       | Print this help message and exit
|General Options
|  -b,--bam FILE REQUIRED          | bam file
|  -g,--genome FILE REQUIRED       | reference genome
|  -a,--anno FILE REQUIRED         | annotation database file
|  -r,--reg FILE                   | valid region to disvover SV
|  -o,--bcfout TEXT=out.bcf        | output bcf file
|  -t,--tsvout TEXT=out.tsv        | output tsv file
|  -s,--svtype INT in [0 - 4]      | SV types to discover,0:INV,1:DEL,2:DUP,3:INS,4:BND
|  -n,--nthread INT in [1 - 20]=8  | number of threads used to process bam
|Threshold Options
|  --min_ref_sep INT=50            | min sv length to compute
|  --max_read_sep INT=10           | max read split mapping pos allowed to compute sv
|  --min_flk_len INT=10            | min flank length needed on each side of breakpoint
|  --min_map_qual INT=1            | min mapping quality of read used to compute sv
|  --min_clip_len INT=20           | min clip length needed to compute sv
|  --min_tra_qual INT=1            | min mapping quality of read pair used to compute sv
|  --min_dup_dp INT=100            | min duplication size needed for an DP to compute sv
|  --min_inv_rpt INT=100           | min inversion size to report
|  --min_del_rpt INT=300           | min deletion size to report
|  --min_dup_rpt INT=100           | min dup size to report

Installation   

1. clone repo   
`git clone https://github.com/vanNul/sver`  

2. compile  
`cd sver`    
`./autogen.sh`   
`./configure --prefix=/path/to/install/dir/`  
`make` 
`make install`

3. execute  
`/path/to/install/dir/sver` 

4. test  
`/path/to/install/dir/sver -b testdata/test.bam -g /Users/wood/Database/hg19/hg19.fa -a testdata/refGene.Anno.sorted.gz -r testdata/valid.bed -o testdata/sv.bcf -t testdata/sv.tsv` 

5. document  
a simple FAQ collection is [here](./doc/README.md)