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
    if(std::abs(pos1 - other.pos1) > 100) return false;
    if(std::abs(pos2 - other.pos2) > 100) return false;
    return true;
}

bool RealnFilter::validCCSeq(const std::string& seq, const std::string& chr1, int32_t& pos1, const std::string& chr2, int32_t& pos2){
    return true;
    std::vector<bam1_t*> alnret;
    mBWA->alignSeq("seq", seq, alnret);
    std::vector<bam1_t*> palnret;
    for(auto& e: alnret){
        if(e->core.flag & (BAM_FSECONDARY | BAM_FUNMAP)){
            bam_destroy1(e);
        }else{
            std::pair<int32_t, int32_t> clip = bamutil::getSoftClipLength(e);
            if(clip.first + clip.second == 0){
                bam_destroy1(e);
            }else{
                palnret.push_back(e);
            }
        }
    }
    if(palnret.size() <= 1){
        for(auto& e: palnret) bam_destroy1(e);
        if(palnret.size()) return true;// one sc
        else return false; // no sc
    }
    // valid scs
    bool valid = false;
    BpPair obp;
    obp.tid1 = bam_name2id(mHeader, chr1.c_str());
    obp.tid2 = bam_name2id(mHeader, chr2.c_str());
    obp.pos1 = pos1;
    obp.pos2 = pos2;
    obp.adjustpt();
    for(uint32_t i = 0; i < palnret.size() - 1; ++i){
        std::vector<BpPair> nbps;
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
        if(lsc){
            BpPair nbp;
            nbp.tid1 = palnret[i]->core.tid;
            nbp.pos1 = palnret[i]->core.pos;
            nbps.push_back(nbp);
        }
        if(rsc){
             BpPair nbp;
             nbp.tid1 = palnret[i]->core.tid;
             nbp.pos1 = r;
             nbps.push_back(nbp);
        }
        std::vector<BpPair> mbps;
        for(uint32_t k = i + 1; k < palnret.size(); ++k){
            data = bam_get_cigar(palnret[k]);
            int rk = palnret[k]->core.pos;
            int lsck = 0;
            int rsck = 0;
            for(uint32_t j = 0; j < palnret[i]->core.n_cigar; ++j){
                uint32_t oplen = bam_cigar_oplen(data[j]);
                int opmask = bam_cigar_op(data[j]);
                switch(bam_cigar_type(opmask)){
                    case 2: case 3:
                        rk += oplen;
                        break;
                    default:
                        break; 
                }
                if(opmask == BAM_CSOFT_CLIP){
                    if(j == 0) lsck = oplen;
                    else rsck = oplen;
                }
            }
            if(lsck){
                BpPair mbp;
                mbp.tid2 = palnret[k]->core.tid;
                mbp.pos2 = palnret[k]->core.pos;
                mbps.push_back(mbp);
            }
            if(rsck){
                BpPair mbp;
                mbp.tid2 = palnret[k]->core.tid;
                mbp.pos2 = rk;
                mbps.push_back(mbp);
            }
            for(uint32_t nc = 0; nc < nbps.size(); ++nc){
                for(uint32_t mc = 0; mc < mbps.size(); ++mc){
                    BpPair sbp;
                    sbp.tid1 = nbps[nc].tid1;
                    sbp.pos1 = nbps[nc].pos1;
                    sbp.tid2 = mbps[mc].tid2;
                    sbp.pos2 = mbps[mc].pos2;
                    sbp.adjustpt();
                    if(obp.agree(sbp)){
                        valid = true;
                        goto gotvalid;
                    }
                }
            }
        }
    }
gotvalid:
    for(auto& e: palnret) bam_destroy1(e);
    return valid;
}
