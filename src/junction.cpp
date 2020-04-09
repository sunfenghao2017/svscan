#include "junction.h"
#include <bamutil.h>

int JunctionMap::insertJunction(const bam1_t* b, bam_hdr_t* h){
    // sa needed
    uint8_t* sa = bam_aux_get(b, "SA");
    if(!sa) return -2;
    // parse non-supplementary alignment record first
    std::vector<Junction> jcvec;
    bool fw = !(b->core.flag & BAM_FREVERSE);
    int32_t refpos = b->core.pos;
    int32_t seqpos = 0, readpos = 0, seqmatch = 0, psc = 0;
    int32_t readStart = refpos;
    uint32_t* cigar = bam_get_cigar(b);
    int32_t seqlen = bamutil::getSeqLen(b);
    int clpst = 0;
    int sccnt = 0;
    for(uint32_t i = 0; i < b->core.n_cigar; ++i){
        int opint = bam_cigar_op(cigar[i]);
        int oplen = bam_cigar_oplen(cigar[i]);
        if(opint == BAM_CMATCH || opint == BAM_CEQUAL || opint == BAM_CDIFF){
            refpos += oplen;
            seqpos += oplen;
        }else if(opint == BAM_CDEL){
            refpos += oplen;
        }else if(opint == BAM_CINS){
            seqpos += oplen;
        }else if(opint == BAM_CSOFT_CLIP){
            ++sccnt;
            int32_t lastSeqPos = seqpos;
            bool scleft = false;
            if(seqpos == 0){
                lastSeqPos += oplen;
                scleft = true;
            }
            seqpos += oplen;
            readpos = (lastSeqPos <= seqlen && !fw) ? seqlen - lastSeqPos : lastSeqPos;
            seqmatch = seqlen - oplen;
            psc = oplen;
            if(oplen > mOpt->filterOpt->minClipLen){
                jcvec.push_back(Junction(fw, scleft, oplen, b->core.tid, readStart, refpos, readpos, seqmatch, b->core.flag & BAM_FREAD1));
            }
        }else if(opint == BAM_CHARD_CLIP){
            return -1;
        }else if(opint == BAM_CREF_SKIP) refpos += oplen;
    }
    if(sccnt > 1) return sccnt;
    clpst = jcvec.size();
    if(clpst != 1) return clpst;
    // parse supplenmentary alignment record then
    if(sa){
        std::string sastr = bam_aux2Z(sa);
        // get optimal SA
        std::vector<std::string> cvs;
        std::vector<std::string> vstr;
        util::split(sastr, cvs, ";");
        if(cvs[1].empty()){
            if(cvs[0].find_first_of("SH") == cvs[0].find_last_of("SH")) sastr = cvs[0];
            else sastr = "";
        }else{
            std::vector<int32_t> mvidx;
            for(uint32_t cvidx = 0; cvidx < cvs.size() - 1; ++cvidx){
                util::split(cvs[cvidx], vstr, ",");
                if(vstr[3].find_first_of("SH") == vstr[3].find_last_of("SH")){
                    mvidx.push_back(cvidx);
                }
            }
            if(mvidx.size() == 1){
                sastr = cvs[mvidx[0]];
            }else{
                sastr = "";
            }
        }
        if(sastr.empty()){
            return false;
        }
        util::split(sastr, vstr, ",");
        int32_t tid = bam_name2id(h, vstr[0].c_str());
        refpos = std::atoi(vstr[1].c_str()) - 1;
        fw = (vstr[2][0] == '+');
        readStart = -1;
        seqpos = 0, readpos = 0, seqmatch = 0;
        std::vector<std::pair<int32_t, char>> pcigar;
        bamutil::parseCigar(vstr[3], pcigar);
        for(auto& e: pcigar){
            char opchr = e.second;
            int oplen = e.first;
            if(opchr == 'M' || opchr == '=' || opchr == 'X'){
                refpos += oplen;
                seqpos += oplen;
            }else if(opchr == 'D'){
                refpos += oplen;
            }else if(opchr == 'I'){
                seqpos += oplen;
            }else if(opchr == 'S'){
                int32_t lastSeqPos = seqpos;
                bool scleft = false;
                if(seqpos == 0){
                    lastSeqPos += oplen;
                    scleft = true;
                }
                seqpos += oplen;
                readpos = (lastSeqPos <= seqlen && !fw) ? seqlen - lastSeqPos : lastSeqPos;
                seqmatch = seqlen - oplen;
                if((seqmatch > mOpt->filterOpt->mMinGoodSRLen || seqmatch <= psc) &&  oplen > mOpt->filterOpt->minClipLen){
                    jcvec.push_back(Junction(fw, scleft, oplen, tid, readStart, refpos, readpos, seqmatch, b->core.flag & BAM_FREAD1));
                }
            }else if(opchr == BAM_CREF_SKIP) refpos += oplen;
        }
    }
    if(jcvec.size() == 2){
        size_t seed = svutil::hashString(bam_get_qname(b));
        auto iter = mJunctionReads.find(seed);
        if(iter == mJunctionReads.end()){
            mJunctionReads[seed] = jcvec;
        }else{
            std::copy(jcvec.begin(), jcvec.end(), std::back_inserter(iter->second));
        }
    }
    return clpst;
}
