#include "dpbamrecord.h"

void DPBamRecord::initClique(int32_t& svStart, int32_t& svEnd, int32_t svt){
    if(svt >= 5){
        int ct = svt - 5;
        if(ct == 0 || ct == 2){// 5to5 and 5to3
            svStart = mCurPos + mCurAlen;
            if(ct == 2) svEnd = mMatePos;
            else svEnd = mMatePos + mMateAlen;
        }else{// 3to3 and 3to5
            svStart = mCurPos;
            if(ct == 3) svEnd = mMatePos + mMateAlen;
            else svEnd = mMatePos;
        }
    }else{
        if(svt == 0){// 5to5 left spanning inversion
            svStart = mMatePos + mMateAlen;
            svEnd = mCurPos + mCurAlen;
        }else if(svt == 1){// 3to3 right spanning inversion
            svStart = mMatePos;
            svEnd = mCurPos;
        }else if(svt == 2){// 5to3 deletion
            svStart = mMatePos + mMateAlen;
            svEnd = mCurPos;
        }else if(svt == 3){// 3to5 duplication
            svStart = mMatePos;
            svEnd = mCurPos + mCurAlen;
        }
    }
}

void DPBamRecord::updateClique(int32_t& svStart, int32_t& svEnd, int32_t svt){
    if(svt >= 5){
        int ct = svt - 5;
        if(ct == 0 || ct == 2){// 5to5 and 5to3
            svStart = std::max(svStart, mCurPos + mCurAlen);
            if(ct == 2){// 5to3
                svEnd = std::min(svEnd, mMatePos);
            }else{// 5to5
                svEnd = std::max(svEnd, mMatePos + mMateAlen);
            }
        }else{// 3to3 and 3to5
            svStart = std::min(svStart, mCurPos);
            if(ct == 3){// 3to5
                svEnd = std::max(svEnd, mMatePos + mMateAlen);
            }else{// 3to3
                svEnd = std::min(svEnd, mMatePos);
            }
        }
    }else{
        if(svt == 0){ // 5to5
            svStart = std::max(svStart, mMatePos + mMateAlen);
            svEnd = std::max(svEnd, mCurPos + mCurAlen);
        }else if(svt == 1){ // 3to3
            svStart = std::min(svStart, mMatePos);
            svEnd = std::min(svEnd, mCurPos);
        }else if(svt == 2){// 5to3
            svStart = std::max(svStart, mMatePos + mMateAlen);
            svEnd = std::min(svEnd, mCurPos);
        }else if(svt == 3){// 3to5
            svStart = std::min(svStart, mMatePos);
            svEnd = std::max(svEnd, mCurPos + mCurAlen);
        }
    }
}

int DPBamRecord::getSVType(const bam1_t* b){
    int svt = -1;
    // check basic SV type
    if(!(b->core.flag & BAM_FREVERSE)){
        if(!(b->core.flag & BAM_FMREVERSE)) svt = 0;
        else{
            if(b->core.pos < b->core.mpos) svt = 2;
            else svt = 3;
        }
    }else{
        if(b->core.flag & BAM_FMREVERSE) svt = 1;
        else{
            if(b->core.pos > b->core.mpos) svt = 2;
            else svt = 3;
        }
    }
    return svt;
}

int DPBamRecord::getSVType(const bam1_t* b, Options* opt){
    int svt = -1;
    // check basic SV type
    if(!(b->core.flag & BAM_FREVERSE)){
        if(!(b->core.flag & BAM_FMREVERSE)) svt = 0;//F...F
        else{
            if(b->core.pos < b->core.mpos) svt = 2; //F...R
            else svt = 3; //R...F
        }
    }else{
        if(b->core.flag & BAM_FMREVERSE) svt = 1;//R...R
        else{
            if(b->core.pos < b->core.mpos) svt = 3;//R...F
            else svt = 2;//F..R
        }
    }
    // refine translocation type
    if(b->core.tid != b->core.mtid){
        if(svt == 0 || svt == 1){
            svt += 5;
            return svt;
        }
        if(b->core.tid < b->core.mtid){
            if(b->core.flag & BAM_FREVERSE) svt = 7;
            else svt = 8;
        }else{
            if(b->core.flag & BAM_FREVERSE) svt = 8;
            else svt = 7;
        }
        return svt;
    }
    // isize checking for DP on same chr
    if(b->core.mpos == b->core.pos) return -1; // NO SV
    if(svt == 2 && std::abs(b->core.isize) < opt->libInfo->mMaxISizeCutoff) return -1; // Too short Del
    if(svt == 3 && std::abs(b->core.pos - b->core.mpos) < opt->filterOpt->mMinDupSize) return -1; // Too short Dup
    return svt;
}

void DPBamRecordSet::cluster(std::vector<DPBamRecord> &dps, SVSet &svs, int32_t svt){
    int32_t origSize = svs.size();
    util::loginfo("Beg clustering DPs for SV type " + std::to_string(svt) + ", all " + std::to_string(dps.size()) + " DPs ");
    std::set<int32_t> clique; // components cluster
    int32_t totdps = dps.size();
    int32_t i = 0, j = 0;
    while(i < totdps){
        clique.clear();
        clique.insert(i);
        j = i + 1;
        while(j < totdps){
            if(dps[j].mCurTid == dps[i].mCurTid &&
               dps[j].mMateTid == dps[i].mMateTid &&
               dps[j].mCurPos - dps[i].mCurPos < mOpt->libInfo->mVarisize &&
               std::abs(dps[j].mMatePos - dps[i].mMatePos) < mOpt->libInfo->mVarisize){
                clique.insert(j);
                ++j;
            }else{
                break;
            }
        }
        searchCliques(clique, dps, svs, svt);
        i = j;
    }
    util::loginfo("End clustering DPs for SV type " + std::to_string(svt) + ", got " + std::to_string(svs.size() - origSize) + " SV candidates.");
}

void DPBamRecordSet::searchCliques(std::set<int32_t>& clique, std::vector<DPBamRecord>& dps, SVSet& svs, int32_t svt){
    auto iter = clique.begin();
    int32_t dpid = *iter;
    int32_t svStart = -1, svEnd = -1, minStart = -1, minEnd = -1, maxStart = -1, maxEnd = -1;
    int32_t chr1 = dps[dpid].mCurTid;
    int32_t chr2 = dps[dpid].mMateTid;
    dps[dpid].initClique(svStart, svEnd, svt);
    minStart = svStart, maxStart = svStart, minEnd = svEnd, maxEnd = svEnd;
    ++iter;
    for(; iter != clique.end(); ++iter){
        dpid = *iter;
        dps[dpid].updateClique(svStart, svEnd, svt);
        minStart = std::min(svStart, minStart);
        maxStart = std::max(svStart, maxStart);
        minEnd = std::min(svEnd, minEnd);
        maxEnd = std::max(svEnd, maxEnd);
    }
    if(clique.size() >= mOpt->filterOpt->mMinSeedDP && validSVSize(svStart, svEnd, svt)){
        SVRecord* svr = new SVRecord();
        svr->mChr1 = chr1;
        svr->mChr2 = chr2;
        svr->mSVStart = svStart;
        svr->mSVEnd = svEnd;
        svr->mPESupport = clique.size();
        int32_t posVar = std::max(std::abs(maxStart - minStart), 50);
        int32_t endVar = std::max(std::abs(maxEnd - minEnd), 50);
        svr->mCiPosLow = -posVar;
        svr->mCiPosHigh = posVar;
        svr->mCiEndLow = -endVar;
        svr->mCiEndHigh = endVar;
        std::vector<uint8_t> mapQV(clique.size(), 0);
        int qIdx = 0;
        for(auto& e : clique) mapQV[qIdx++] = dps[e].mMapQual;
        svr->mPEMapQuality = statutil::median(mapQV);
        svr->mSRSupport = 0;
        svr->mSRAlignQuality = 0;
        svr->mPrecise = 0;
        svr->mSVT = svt;
        svr->mAlnInsLen = 0;
        svr->mHomLen = 0;
        svs.push_back(svr);
        int32_t svid = svs.size();
        for(auto& e: clique) dps[e].mSVID = svid;
    }
}

void DPBamRecordSet::cluster(SVSet& svs){
    for(auto& e : mOpt->SVTSet){
        if(!mDPs[e].empty()){
            std::sort(mDPs[e].begin(), mDPs[e].end());
            cluster(mDPs[e], svs, e);
        }
    }
}
