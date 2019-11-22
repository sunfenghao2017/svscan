#include "realnfilter.h"

void BpPair::adjustpt(){
    if(tid1 > tid2){
        std::swap(tid1, tid2);
        std::swap(pos1, pos2);
    }
    if(tid1 == tid2 && pos1 > pos2){
        std::swap(pos1, pos2);
    }
}

bool BpPair::agree(const BpPair& other){
    if(tid1 != other.tid1) return false;
    if(tid2 != other.tid2) return false;
    if(std::abs(pos1 - other.pos1) > 10) return false;
    if(std::abs(pos2 - other.pos2) > 10) return false;
    return true;
}

bool RealnFilter::validCCSeq(const std::string& seq, const std::string& chr1, int32_t pos1, const std::string& chr2, int32_t pos2){
    std::vector<bam1_t*> alnret;
    mBWA->alignSeq("seq", seq, alnret);
    std::vector<bam1_t*> palnret;
    for(auto& e: alnret){
        if((e->core.flag & (BAM_FSECONDARY | BAM_FUNMAP)) || e->core.qual <= 0){
            bam_destroy1(e);
        }else{
            std::pair<int32_t, int32_t> clip = bamutil::getSoftClipLength(e);
            if(clip.first + clip.second == 0){// no clip
                bam_destroy1(e);
            }else if(clip.first ^ clip.second){// good, one clip
                palnret.push_back(e);
            }else if(clip.first && clip.second){
                if((clip.first < 10 && clip.second > 10) || (clip.first > 10 && clip.second < 10)){
                    palnret.push_back(e);
                }else{
                    bam_destroy1(e);
                }
            }
        }
    }
    if(palnret.size() != 2){
        for(auto& e: palnret) bam_destroy1(e);
        if(palnret.size()) return true;
        else return false;
    }
    bool valid = true;
    BpPair bp;
    for(uint32_t i = 0; i < 2; ++i){
        uint32_t* data = bam_get_cigar(palnret[i]);
        int r = palnret[i]->core.pos;
        int lsc = 0;
        int rsc = 0;
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
                else rsc = oplen;
            }
        }
        if(i){
            bp.tid2 = palnret[i]->core.tid;
            if(lsc < 10) bp.pos2 = palnret[i]->core.pos;
            else bp.pos2 = r;
        }else{
            bp.tid1 = palnret[i]->core.tid;
            if(lsc < 10) bp.pos1 = palnret[i]->core.pos;
            else bp.pos1 = r;
        }
    }
    for(auto& e: palnret) bam_destroy1(e);
    if(!valid) return false;
    BpPair obp;
    obp.tid1 = bam_name2id(mHeader, chr1.c_str());
    obp.tid2 = bam_name2id(mHeader, chr2.c_str());
    obp.pos1 = pos1;
    obp.pos2 = pos2;
    bp.adjustpt();
    obp.adjustpt();
    return bp.agree(obp);
}
