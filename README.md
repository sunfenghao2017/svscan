program: sver
version: 0.0.0
updated: 21:21:42
Usage: sver [OPTIONS]

Options:
  -h,--help                         Print this help message and exit

General:
  -b,--bam FILE REQUIRED            bam file
  -g,--genome FILE REQUIRED         reference genome
  -a,--anno FILE REQUIRED           annotation database file
  -r,--reg FILE                     valid region to disvover SV
  -o,--bcfout TEXT=out.bcf          output bcf file
  -t,--tsvout TEXT=out.tsv          output tsv file
  -s,--svtype INT in [0 - 4]        SV types to discover,0:INV,1:DEL,2:DUP,3:INS,4:BND
  -n,--nthread INT in [1 - 20]=8    number of threads used to process bam
