#include "svbam.h"

void SVBAMOpt::getBam(){
    std::string tobam = obam + ".tmp.unsorted.bam";
    std::string tobam2 = obam2 + ".tmp.unsorted.bam";
    const char* idt = "ZF";
    samFile* ifp = sam_open(ibam.c_str(), "r");
    samFile* ofp = sam_open(tobam.c_str(), "w");
    samFile* ofp2 = NULL;
    bam_hdr_t* h = sam_hdr_read(ifp);
    assert(sam_hdr_write(ofp, h) >= 0);
    if(outalt){
        ofp2 = sam_open(tobam2.c_str(), "w");
        assert(sam_hdr_write(ofp2, h) >= 0);
    }
    bam1_t* b = bam_init1();
    std::vector<bam1_t*> realn;
    while(sam_read1(ifp, h, b) >= 0){
        if(bam_aux2i(bam_aux_get(b, idt)) == svid){
            assert(sam_write1(ofp, h, b) >= 0);
        }
        if(outalt){
            UnalignedSeq us;
            us.mName = bamutil::getQName(b);
            us.mSeq = bamutil::getSeq(b);
            realn.clear();
            obwa->alignSeq(us, realn);
            for(auto& e: realn){
                assert(sam_write1(ofp2, h, e) >= 0);
                bam_destroy1(e);
            }
        }
    }
    sam_close(ifp);
    sam_close(ofp);
    if(outalt) sam_close(ofp2);
    bam_hdr_destroy(h);
    bam_destroy1(b);
    std::string sortCMD = "samtools sort -o " + obam + " " + tobam;
    system(sortCMD.c_str());
    remove(tobam.c_str());
    if(outalt){
        sortCMD = "samtools sort -o " + obam2 + " " + tobam2;
        system(sortCMD.c_str());
        remove(tobam2.c_str());
        assert(sam_index_build(obam2.c_str(), 14) == 0);
    }
    assert(sam_index_build(obam.c_str(), 14) == 0);
}
