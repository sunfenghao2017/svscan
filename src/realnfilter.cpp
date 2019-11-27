#include "realnfilter.h"

void BpPair::adjustpt(){
    if(tid1 > tid2){
        std::swap(tid1, tid2);
        std::swap(pos1, pos2);
        swapped = true;
    }
    if(tid1 == tid2 && pos1 > pos2){
        std::swap(pos1, pos2);
        swapped = true;
    }
}

bool BpPair::agree(const BpPair& other){
    if(tid1 != other.tid1) return false;
    if(tid2 != other.tid2) return false;
    return true;
}

bool RealnFilter::validCCSeq(const std::string& seq, const std::string& chr1, int32_t& pos1, const std::string& chr2, int32_t& pos2){
    std::vector<bam1_t*> alnret;
    mBWA->alignSeq("seq", seq, alnret);
    std::vector<bam1_t*> palnret;
    for(auto& e: alnret){
        if(e->core.flag & (BAM_FSECONDARY | BAM_FUNMAP)){
            bam_destroy1(e);
        }else{
            std::pair<int32_t, int32_t> clip = bamutil::getSoftClipLength(e);
            if((clip.first > 0) ^ (clip.second > 0)) palnret.push_back(e);
            else bam_destroy1(e);
        }
    }
    if(palnret.size() != 2){
        for(auto& e: palnret) bam_destroy1(e);
        if(palnret.size() == 1) return true;// one sc
        else return false; // more than two sc
    }
    // valid scs
    BpPair obp;
    obp.tid1 = bam_name2id(mHeader, chr1.c_str());
    obp.tid2 = bam_name2id(mHeader, chr2.c_str());
    obp.pos1 = pos1;
    obp.pos2 = pos2;
    obp.adjustpt();
    BpPair nbp;
    for(uint32_t i = 0; i < 2; ++i){
        std::vector<BpPair> nbps;
        uint32_t* data = bam_get_cigar(palnret[i]);
        int r = palnret[i]->core.pos;
        int lsc = 0;
        for(uint32_t j = 0; j < palnret[i]->core.n_cigar; ++j){
            uint32_t oplen = bam_cigar_oplen(data[j]);
            int opmask = bam_cigar_op(data[j]);
            switch(bam_cigar_type(opmask)){
                case 2: case 3:
                    r += oplen;
                    break;
                default:
                    break; 
            }
            if(opmask == BAM_CSOFT_CLIP){
                if(j == 0) lsc = oplen;
            }
        }
        if(i){
            nbp.tid1 = palnret[i]->core.tid;
            if(lsc) nbp.pos1 = palnret[i]->core.pos;
            else nbp.pos1 = r;
        }else{
            nbp.tid2 = palnret[i]->core.tid;
            if(lsc) nbp.pos2 = palnret[i]->core.pos;
            else nbp.pos2 = r;
        }
    }
    nbp.adjustpt();
    for(auto& e: palnret) bam_destroy1(e);
    bool valid = obp.agree(nbp);
    if(valid){
        if(obp.swapped){
            pos2 = std::min(nbp.pos1, (int32_t)mHeader->target_len[nbp.tid1]);
            pos1 = std::min(nbp.pos2, (int32_t)mHeader->target_len[nbp.tid2]);
        }else{
            pos2 = std::min(nbp.pos2, (int32_t)mHeader->target_len[nbp.tid2]);
            pos1 = std::min(nbp.pos1, (int32_t)mHeader->target_len[nbp.tid1]);
        }
    }
    return valid;
}
