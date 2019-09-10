#include <htslib/sam.h>
#include <iostream>
#include <cassert>
#include <string>
#include <vector>

int main(int argc, char** argv){
    if(argc < 3){
        std::cout << argv[0] << " <sv.bam> <idv> [idv.bam] " << std::endl;
        return 0;
    }
    
    char* ibam = argv[1];
    int idv = std::atoi(argv[2]);
    std::string obam = std::string("sv.") + argv[2] + std::string(".bam");
    if(argc > 3) obam = argv[3];
    std::string tobam = obam + ".tmp.unsorted.bam";

    const char* idt = "ZF";
    samFile* ifp = sam_open(ibam, "r");
    samFile* ofp = sam_open(tobam.c_str(), "w");
    bam_hdr_t* h = sam_hdr_read(ifp);
    assert(sam_hdr_write(ofp, h) >= 0);
    bam1_t* b = bam_init1();
    while(sam_read1(ifp, h, b) >= 0){
        if(bam_aux2i(bam_aux_get(b, idt)) == idv){
            assert(sam_write1(ofp, h, b) >= 0);
        }
    }
    sam_close(ifp);
    sam_close(ofp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
    std::string sortCMD = "samtools sort -o " + obam + " " + tobam;
    system(sortCMD.c_str());
    remove(tobam.c_str());
    assert(sam_index_build(obam.c_str(), 14) == 0);
}
