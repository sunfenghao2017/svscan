#include "dpbamrecord.h"

void DPBamRecord::initClique(int32_t& svStart, int32_t& svEnd, int32_t& wiggle, Options* opt, int32_t svt){
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
    wiggle = opt->libInfo->mMaxNormalISize;
}

bool DPBamRecord::updateClique(int32_t& svStart, int32_t& svEnd, int32_t& wiggle, int32_t svt){
    if(svt >= 5){
        int ct = svt - 5;
        int32_t newSVStart;
        int32_t newSVEnd;
        int32_t newWiggle = wiggle;
        if(ct == 0 || ct == 2){// 5to5 and 5to3
            newSVStart = std::max(svStart, mCurPos + mCurAlen);
            newWiggle -= (newSVStart - svStart);
            if(ct == 2){// 5to3
                newSVEnd = std::min(svEnd, mMatePos);
                newWiggle -= (svEnd - newSVEnd);
            }else{// 5to5
                newSVEnd = std::max(svEnd, mMatePos + mMateAlen);
                newWiggle -= (newSVEnd - svEnd);
            }
        }else{// 3to3 and 3to5
            newSVStart = std::min(svStart, mCurPos);
            newWiggle -= (svStart - newSVStart);
            if(ct == 3){// 3to5
                newSVEnd = std::max(svEnd, mMatePos + mMateAlen);
                newWiggle -= (newSVEnd - svEnd);
            }else{// 3to3
                newSVEnd = std::min(svEnd, mMatePos);
                newWiggle -= (svEnd - newSVEnd);
            }
        }
        // check whether this still a valid translocation cluster
        if(newWiggle > 0){
            svStart = newSVStart;
            svEnd = newSVEnd;
            wiggle = newWiggle;
            return true;
        }
        return false;
    }else{
        int32_t newSVStart = svStart;
        int32_t newSVEnd = svEnd;
        int32_t newWiggle = wiggle;
        if(svt == 0 || svt == 1){// inversion
            if(svt == 0){// 5to5
                newSVStart = std::max(svStart, mMatePos + mMateAlen);
                newSVEnd = std::max(svEnd, mCurPos + mCurAlen);
                newWiggle -= std::max(newSVStart - svStart, newSVEnd - svEnd);
            }else{// 3to3
                newSVStart = std::min(svStart, mMatePos);
                newSVEnd = std::min(svEnd, mCurPos);
                newWiggle -=  std::max(svStart - newSVStart, svEnd - newSVEnd);
            }
        }else if(svt == 2){// deletion
            newSVStart = std::max(svStart, mMatePos + mMateAlen);
            newSVEnd = std::min(svEnd, mCurPos);
            newWiggle -= std::max(newSVStart - svStart,  svEnd - newSVEnd);
        }else if(svt == 3){// duplication
            newSVStart = std::min(svStart, mMatePos);
            newSVEnd = std::max(svEnd, mCurPos + mCurAlen);
            newWiggle -= std::max(newSVEnd - svEnd, svStart - newSVStart);
        }
        // check whether new inversion size agree with all pairs
        if(newSVStart < newSVEnd && newWiggle >= 0){
            svStart = newSVStart;
            svEnd = newSVEnd;
            wiggle = newWiggle;
            return true;
        }
        return false;
    }
    return false;
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
    if(svt == 2 && std::abs(b->core.isize) < opt->filterOpt->mMinInDelSize) return -1; // Too short Del
    if(svt == 3 && std::abs(b->core.pos - b->core.mpos) < opt->filterOpt->mMinDupSize) return -1; // Too short Dup
    return svt;
}

void DPBamRecordSet::cluster(std::vector<DPBamRecord> &dps, SVSet &svs, int32_t svt){
    if(dps.empty()) return;
    // Sort DPBamRecords
    std::sort(dps.begin(), dps.end());
    // Components
    std::vector<int32_t> comp = std::vector<int32_t>(dps.size(), 0);
    // Edge lists for each component
    std::map<int32_t, std::vector<EdgeRecord>> compEdge;
    int32_t compNum = 0;
    size_t lastConnectedNodesEnd = 0;
    size_t lastConnectedNodesBeg = 0;
    for(uint32_t i = 0; i < dps.size(); ++i){
        // Safe to clean the graph
        if(i > lastConnectedNodesEnd){
            // Clean edge lists
            if(!compEdge.empty()){
                searchCliques(compEdge, dps, svs, svt);
                compEdge.clear();
            }
        }
        // Try to expand components
        int32_t mincrd = dps[i].minCoord();
        int32_t maxcrd = dps[i].maxCoord();
        for(uint32_t j = i + 1; j < dps.size(); ++j){
            // Check that mate chr agree
            if(dps[i].mCurTid != dps[j].mCurTid) continue;
            // Check two DP leftmost mapping position falling in reasonable range
            if(std::abs(dps[j].minCoord() + dps[j].mCurAlen - mincrd) > mOpt->libInfo->mVarisize) break;
            // Check two DP rightmost mapping position falling in reasonable range
            if(std::abs(dps[j].maxCoord() - maxcrd) + mOpt->libInfo->mReadLen > mOpt->libInfo->mMaxNormalISize) break;
            // Update last connected node
            if(j > lastConnectedNodesEnd) lastConnectedNodesEnd = j;
            // Assign components
            int32_t compIndex = 0;
            if(!comp[i]){
                if(!comp[j]){// Neither vertex has component assigned
                    compIndex = ++compNum;
                    comp[i] = compIndex;
                    comp[j] = compIndex;
                    compEdge.insert(std::make_pair(compIndex, std::vector<EdgeRecord>()));
                }else{// Only one vertex has component assigned
                    compIndex = comp[j];
                    comp[i] = compIndex;
                }
            }else{
                if(!comp[j]){// Only one veretx has comopnent assigned
                    compIndex = comp[i];
                    comp[j] = compIndex;
                }else{// Both vertices have component assigned
                    if(comp[i] == comp[j]) compIndex = comp[i];
                    else{// Merge components
                        compIndex = comp[i];
                        int32_t otherIndex = comp[j];
                        if(otherIndex < compIndex){
                            compIndex = comp[j];
                            otherIndex = comp[i];
                        }
                        // Re-label other index
                        for(size_t k = lastConnectedNodesBeg; k <= lastConnectedNodesEnd; ++k){
                            if(otherIndex == comp[i]) comp[i] = compIndex;
                        }
                        // Merge edge lists
                        auto compIdxIter = compEdge.find(compIndex);
                        auto otherIdxIter = compEdge.find(otherIndex);
                        compIdxIter->second.insert(compIdxIter->second.end(), otherIdxIter->second.begin(), otherIdxIter->second.end());
                        compEdge.erase(otherIdxIter);
                    }
                }
            }
            // Append new edge
            auto compEdgeIter = compEdge.find(compIndex);
            if(compEdgeIter->second.size() < mOpt->filterOpt->mGraphPruning){
                int32_t weight = std::abs(std::abs(dps[i].minCoord() - dps[j].minCoord()) - std::abs(dps[i].maxCoord() - dps[j].maxCoord()));
                compEdgeIter->second.push_back(EdgeRecord(i, j, weight));
            }
        }
    }
    if(!compEdge.empty()){
        searchCliques(compEdge, dps, svs, svt);
        compEdge.clear();
    }
}

void DPBamRecordSet::searchCliques(std::map<int32_t, std::vector<EdgeRecord>>& compEdge, std::vector<DPBamRecord>& dps, SVSet& svs, int32_t svt){
    // Iterate all components
    for(auto compIter = compEdge.begin(); compIter != compEdge.end(); ++compIter){
        // Sort edges by weight
        std::sort(compIter->second.begin(), compIter->second.end());
        // Find a large clique
        auto edgeIter = compIter->second.begin();
        std::set<int32_t> clique, incompatible;
        int32_t svStart = -1;
        int32_t svEnd = -1;
        int32_t wiggle = 0;
        int32_t chr1 = dps[edgeIter->mSource].mCurTid;
        int32_t chr2 = dps[edgeIter->mSource].mMateTid;
        dps[edgeIter->mSource].initClique(svStart, svEnd, wiggle, mOpt, svt);
        if(chr1 == chr2 && svStart >= svEnd) continue; // do not support SV
        clique.insert(edgeIter->mSource);
        // Grow the clique
        ++edgeIter;
        for(; edgeIter != compIter->second.end(); ++edgeIter){
            int32_t v;
            if(clique.find(edgeIter->mSource) == clique.end() && clique.find(edgeIter->mTarget) != clique.end()){
                v = edgeIter->mSource;
            }else if(clique.find(edgeIter->mSource) != clique.end() && clique.find(edgeIter->mTarget) == clique.end()){
                v = edgeIter->mTarget;
            }else continue;
            if(incompatible.find(v) != incompatible.end()) continue;
            if(dps[v].updateClique(svStart, svEnd, wiggle, svt)) clique.insert(v);
            else incompatible.insert(v);
        }
        if(clique.size() > 1 && validSVSize(svStart, svEnd, svt)){
            SVRecord svr;
            svr.mChr1 = chr1;
            svr.mChr2 = chr2;
            svr.mSVStart = svStart + 1;
            svr.mSVEnd = svEnd + 1;
            svr.mPESupport = clique.size();
            int32_t ciwiggle = std::max(std::abs(wiggle), 50);
            svr.mCiPosLow = -ciwiggle;
            svr.mCiPosHigh = ciwiggle;
            svr.mCiEndLow = -ciwiggle;
            svr.mCiEndHigh = ciwiggle;
            std::vector<uint8_t> mapQV(clique.size(), 0);
            int qIdx = 0;
            for(auto& e : clique) mapQV[qIdx++] = dps[e].mMapQual;
            svr.mPEMapQuality = statutil::median(mapQV);
            svr.mSRSupport = 0;
            svr.mSRAlignQuality = 0;
            svr.mPrecise = 0;
            svr.mSVT = svt;
            svr.mInsLen = 0;
            svr.mHomLen = 0;
            svs.push_back(svr);
        }
    }
}

void DPBamRecordSet::cluster(SVSet& svs){
    for(auto& e : mOpt->SVTSet){
        cluster(mDPs[e], svs, e);
    }
}
