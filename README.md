program: sver  
version: 0.0.0  
updated: 22:01:59 Aug  6 2019  
Usage: sver [OPTIONS]  

|Options                            | Explanations
|-----------------------------------|--------------------------------------------------
|  -h,--help                        | Print this help message and exit
|General:                           
|  -b,--bam FILE REQUIRED           | bam file
|  -g,--genome FILE REQUIRED        | reference genome
|  -a,--anno FILE REQUIRED          | annotation database file
|  -r,--reg FILE                    | valid region to disvover SV
|  -o,--bcfout TEXT=out.bcf         | output bcf file
|  -t,--tsvout TEXT=out.tsv         | output tsv file
|  -s,--svtype INT in [0 - 4]       | SV types to discover,0:INV,1:DEL,2:DUP,3:INS,4:BND
|  -n,--nthread INT in [1 - 20]=8   | number of threads used to process bam

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
