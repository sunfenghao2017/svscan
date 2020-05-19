#include "svbam.h"

void SVBAMOpt::getBam(){
    std::vector<bam1_t*> brso, brso2;
    const char* idt = "ZF";
    samFile* ifp = sam_open(ibam.c_str(), "r");
    samFile* ofp = sam_open(obam.c_str(), "w");
    samFile* ofp2 = NULL;
    bam_hdr_t* h = sam_hdr_read(ifp);
    bam_hdr_t* h2 = NULL;
    assert(sam_hdr_write(ofp, h) >= 0);
    if(outalt){
        ofp2 = sam_open(obam2.c_str(), "w");
        h2 = obwa->getBamHeader();
        assert(sam_hdr_write(ofp2, h2) >= 0);
    }
    bam1_t* b = bam_init1();
    std::vector<bam1_t*> realn;
    while(sam_read1(ifp, h, b) >= 0){
        if(bam_aux2i(bam_aux_get(b, idt)) == svid){
            brso.push_back(b);
            if(outalt){
                UnalignedSeq us;
                us.mName = bamutil::getQName(b);
                us.mSeq = bamutil::getSeq(b);
                realn.clear();
                obwa->alignSeq(us, realn);
                for(auto& e: realn){
                    if(!(e->core.flag& BAM_FSECONDARY)) brso2.push_back(e);
                    else bam_destroy1(e);
                }
            }
            b = bam_init1();
        }
    }
    sam_close(ifp);
    bam_destroy1(b);
    // sort
    std::sort(brso.begin(), brso.end(), BamComp());
    for(auto& e: brso){
        assert(sam_write1(ofp, h, e) >= 0);
        bam_destroy1(e);
    }
    sam_close(ofp);
    assert(sam_index_build(obam.c_str(), 0) == 0);
    if(outalt){
        std::sort(brso2.begin(), brso2.end(), BamComp());
        for(auto& e: brso2){
            assert(sam_write1(ofp2, h, e) >= 0);
            bam_destroy1(e);
        }
        sam_close(ofp2);
        assert(sam_index_build(obam2.c_str(), 0) == 0);
    }
    sam_hdr_destroy(h);
}
