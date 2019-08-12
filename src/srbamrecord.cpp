#include "srbamrecord.h"

void SRBamRecordSet::classifyJunctions(JunctionMap* jctMap){
    util::loginfo("Start classifing SRs into various SV candidates");
    if(!jctMap->mSorted) jctMap->sortJunctions();
    int svtIdx = 0;
    int32_t rst = -1;
    bool hasSA = false;
    for(auto iter = jctMap->mJunctionReads.begin(); iter != jctMap->mJunctionReads.end(); ++iter){
        // find insertion candidates which have no supplementary alignments
        hasSA = false;
        for(uint32_t i = 0; i < iter->second.size(); ++i){
            if(iter->second[i].mRstart > 0){
                hasSA = true;
                break;
            }
        }
        if(!hasSA){
            for(uint32_t i = 0; i < iter->second.size(); ++i){
                svtIdx = 4;
                mSRs[svtIdx].push_back(SRBamRecord(iter->second[i].mRefidx,
                                                   iter->second[i].mRefpos,
                                                   iter->second[i].mRefidx,
                                                   iter->second[i].mRefpos,
                                                   rst,
                                                   std::abs(iter->second[i].mSeqpos - iter->second[i].mSeqpos),
                                                   iter->first));

            }
            continue;
        }
        for(uint32_t i = 0; i < iter->second.size(); ++i){
            for(uint32_t j = i + 1; j < iter->second.size(); ++j){
                if(iter->second[i].mRstart > 0 && iter->second[j].mRstart > 0) continue;
                if(iter->second[i].mSCLen && iter->second[j].mSCLen && 
                   (iter->second[j].mSCLen < iter->second[i].mSeqmatch - mOpt->filterOpt->mMaxReadSep ||
                   iter->second[i].mSCLen < iter->second[j].mSeqmatch - mOpt->filterOpt->mMaxReadSep)){
                    continue;
                }
                // get read starting mapping position
                rst = iter->second[i].mRstart;
                if(rst == -1) rst = iter->second[j].mRstart;
                // check possible translocation split read
                if(iter->second[j].mRefidx != iter->second[i].mRefidx){
                    int32_t littleChrIdx = i;
                    int32_t largerChrIdx = j;
                    if(iter->second[j].mRefidx < iter->second[i].mRefidx){
                        littleChrIdx = j;
                        largerChrIdx = i;
                    }
                    svtIdx = 0;
                    if(iter->second[littleChrIdx].mForward == iter->second[largerChrIdx].mForward){
                        // Same direction, opposing soft-clips
                        if(iter->second[littleChrIdx].mSCleft != iter->second[largerChrIdx].mSCleft){
                            if(iter->second[littleChrIdx].mSCleft){
                                svtIdx = 7; // littleChr on 3' part
                            }else{
                                svtIdx = 8; // littleChr on 5' part
                            }
                        }
                    }else{
                        //opposing direction, same doft-clips
                        if(iter->second[littleChrIdx].mSCleft == iter->second[largerChrIdx].mSCleft){
                            if(iter->second[littleChrIdx].mSCleft){
                                svtIdx = 6; // 3to3 connection
                            }else{
                                svtIdx = 5; // 5to5 connection
                            }
                        }
                    }
                    if(svtIdx && mOpt->SVTSet.find(svtIdx) != mOpt->SVTSet.end()){
                        mSRs[svtIdx].push_back(SRBamRecord(iter->second[largerChrIdx].mRefidx,
                                                        iter->second[largerChrIdx].mRefpos,
                                                        iter->second[littleChrIdx].mRefidx,
                                                        iter->second[littleChrIdx].mRefpos,
                                                        rst,
                                                        std::abs(iter->second[j].mSeqpos - iter->second[i].mSeqpos),
                                                        iter->first));
                    }
                }else{
                    svtIdx = -1;
                    int32_t leftPart = i;
                    int32_t rightPart = j;
                    if(iter->second[j].mRefpos <= iter->second[i].mRefpos){
                        leftPart = j;
                        rightPart = i;
                    }
                    if(iter->second[j].mForward == iter->second[i].mForward && // same direction
                       iter->second[j].mSCleft != iter->second[i].mSCleft  &&  // opposing soft-clips
                       std::abs(iter->second[j].mRefpos - iter->second[i].mRefpos) < mOpt->filterOpt->mMaxReadSep){// breakpoint close
                        svtIdx = 4; 
                    }else if(iter->second[j].mForward == iter->second[i].mForward && // same direction
                             iter->second[j].mSCleft != iter->second[i].mSCleft && // opposing soft-clips
                             std::abs(iter->second[j].mRefpos - iter->second[i].mRefpos) >= mOpt->filterOpt->mMinRefSep){// breakpoint faraway
                        if(iter->second[leftPart].mSCleft) svtIdx = 3; // left part leading soft-clip, duplication
                        else svtIdx = 2; // left part tailing soft-clip, deletion
                    }else if(iter->second[j].mForward != iter->second[i].mForward && // opposing direction
                             iter->second[j].mSCleft == iter->second[i].mSCleft && // same soft-clips
                             std::abs(iter->second[j].mRefpos - iter->second[i].mRefpos) >= mOpt->filterOpt->mMinRefSep){// breakpoint farway
                        if(iter->second[j].mSCleft) svtIdx = 1; // 3to3 right spanning inversion breakpoint
                        else svtIdx = 0; // 5to5 left spanning inversion breakpoint
                    }
                    if(svtIdx != -1 && mOpt->SVTSet.find(svtIdx) != mOpt->SVTSet.end()){
                        mSRs[svtIdx].push_back(SRBamRecord(iter->second[leftPart].mRefidx,
                                                           iter->second[leftPart].mRefpos,
                                                           iter->second[rightPart].mRefidx,
                                                           iter->second[rightPart].mRefpos,
                                                           rst,
                                                           std::abs(iter->second[j].mSeqpos - iter->second[i].mSeqpos),
                                                           iter->first));
                    }
                }
            }
        }
    }
    util::loginfo("Finish classifing SRs into various SV candidates");
}

void SRBamRecordSet::cluster(std::vector<SRBamRecord>& srs, SVSet& svs, int32_t svt){
    util::loginfo("Start clustering SRs for SV type:" + std::to_string(svt));
    for(auto&refIdx : mOpt->svRefID){
        // Components assigned marker
        std::vector<int32_t> comp = std::vector<int32_t>(srs.size(), 0);
        int32_t compNum = 0;
        Cluster compEdge; // component clusters
        // Construct graphs
        size_t lastConnectedNodesEnd = 0;
        size_t lastConnectedNodesBeg = 0;
        for(uint32_t i = 0; i < srs.size(); ++i){
            if(srs[i].mChr1 != refIdx) continue;
            // Safe to clean the graph ?
            if(i > lastConnectedNodesEnd){
                // Clean edge lists
                if(!compEdge.empty()){
                    searchCliques(compEdge, srs, svs, svt);
                    lastConnectedNodesBeg = lastConnectedNodesEnd;
                    compEdge.clear();
                }
            }
            // Search possible connectable node
            for(uint32_t j = i + 1; j < srs.size(); ++j){
                if(srs[j].mChr1 != refIdx) continue; // same chr
                if(srs[j].mPos1 - srs[i].mPos1 > mOpt->filterOpt->mMaxReadSep) break; // mapping position in valid range
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
                    if(!comp[j]){// Only one vertex has component assigned
                        compIndex = comp[i];
                        comp[j] = compIndex;
                    }else{// Both vertices have components assigned, then merge these components
                        if(comp[i] == comp[j]) compIndex = comp[i];
                        else{// Merge components
                            compIndex = comp[i];
                            int32_t otherIndex = comp[j];
                            if(otherIndex < compIndex){
                                compIndex = comp[j];
                                otherIndex = comp[i];
                            }
                            // Re-label other index
                            for(uint32_t k = lastConnectedNodesBeg; k <= lastConnectedNodesEnd; ++k){
                                if(otherIndex == comp[k]){
                                    comp[k] = compIndex;
                                }
                            }
                            // Merge edge list
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
                    // Breakpoint distance
                    int32_t weight = std::abs(srs[j].mPos2 - srs[i].mPos2) + std::abs(srs[j].mPos1 - srs[i].mPos1);
                    compEdgeIter->second.push_back(EdgeRecord(i, j, weight));
                }
            }
        }
        // Search cliques
        if(!compEdge.empty()){
            searchCliques(compEdge, srs, svs, svt);
            compEdge.clear();
        }
    }
    util::loginfo("Finish clustering SRs for SV type:" + std::to_string(svt));
}

void SRBamRecordSet::searchCliques(Cluster& compEdge, std::vector<SRBamRecord>& srs, SVSet& svs, int32_t svt){
    // Iterate all components
    for(auto compIter = compEdge.begin(); compIter != compEdge.end(); ++compIter){
        // Sort edges by weight
        std::sort(compIter->second.begin(), compIter->second.end());
        auto edgeIter = compIter->second.begin();
        // Find a large clique
        std::set<int32_t> clique, incompatible;
        // Initialization clique
        clique.insert(edgeIter->mSource);
        int32_t chr1 = srs[edgeIter->mSource].mChr1;
        int32_t chr2 = srs[edgeIter->mSource].mChr2;
        int32_t ciposlow = srs[edgeIter->mSource].mPos1;
        uint64_t pos1 = srs[edgeIter->mSource].mPos1;
        int32_t ciposhigh = srs[edgeIter->mSource].mPos1;
        int32_t ciendlow = srs[edgeIter->mSource].mPos2;
        uint64_t pos2 = srs[edgeIter->mSource].mPos2;
        int32_t ciendhigh = srs[edgeIter->mSource].mPos2;
        int32_t inslen = srs[edgeIter->mSource].mInslen;
        ++edgeIter;
        // Grow clique
        for(; edgeIter != compIter->second.end(); ++edgeIter){
            // Find next best edge for extension
            int32_t v;
            if(clique.find(edgeIter->mSource) == clique.end() && clique.find(edgeIter->mTarget) != clique.end()){
                v = edgeIter->mSource;
            }else if(clique.find(edgeIter->mSource) != clique.end() && clique.find(edgeIter->mTarget) == clique.end()){
                v = edgeIter->mTarget;
            }else continue;
            if(incompatible.find(v) != incompatible.end()) continue;
            // Try to update clique with this vertex
            int32_t newCiPosLow = std::min(srs[v].mPos1, ciposlow);
            int32_t newCiPosHigh = std::max(srs[v].mPos1, ciposhigh);
            int32_t newCiEndLow = std::min(srs[v].mPos2, ciendlow);
            int32_t newCiEndHigh = std::max(srs[v].mPos2, ciendhigh);
            if((newCiPosHigh - newCiPosLow) < mOpt->filterOpt->mMaxReadSep &&
               (newCiEndHigh - newCiEndLow) < mOpt->filterOpt->mMaxReadSep){// Accept new vertex
                clique.insert(v);
                ciposlow = newCiPosLow;
                pos1 += srs[v].mPos1;
                ciposhigh = newCiPosHigh;
                ciendlow = newCiEndLow;
                pos2 += srs[v].mPos2;
                ciendhigh = newCiEndHigh;
                inslen += srs[v].mInslen;
              }else incompatible.insert(v);
        }
        // At least 2 split read support
        if(clique.size() > 1){
            int32_t svStart = pos1/clique.size();
            int32_t svEnd = pos2/clique.size();
            int32_t svISize = inslen/clique.size();
            int32_t svid = svs.size();
            SVRecord svr;
            svr.mChr1 = chr1;
            svr.mSVStart = svStart;
            svr.mChr2 = chr2;
            svr.mSVEnd = svEnd;
            svr.mCiPosLow = ciposlow - svStart;
            svr.mCiPosHigh = ciposhigh - svStart;
            svr.mCiEndLow = ciendlow - svEnd;
            svr.mCiEndHigh = ciendhigh - svEnd;
            svr.mAlnInsLen = svISize;
            svr.mID = svid;
            svr.mSVT = svt;
            svs.push_back(svr);
            // Reads assigned
            for(auto& e : clique) srs[e].mSVID = svid;
        }
    }
}

void SRBamRecordSet::cluster(SVSet& svs){
    if(!mSorted) sortSRs();
    for(auto& svt : mOpt->SVTSet){
        cluster(mSRs[svt], svs, svt);
    }
}

void SRBamRecordSet::assembleSplitReads(SVSet& svs){
    // Open file handles
    samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
    hts_set_fai_filename(fp, mOpt->genome.c_str());
    hts_idx_t* idx = sam_index_load(fp, mOpt->bamfile.c_str());
    bam_hdr_t* hdr = sam_hdr_read(fp);
    faidx_t* fai = fai_load(mOpt->genome.c_str());
    bam1_t* b = bam_init1();
    // Construct bool filter of SR mapping positions
    for(uint32_t svt = 0; svt < mSRs.size(); ++svt){
        for(uint32_t i = 0; i < mSRs[svt].size(); ++i){
            if(mSRs[svt][i].mSVID == -1 || mSRs[svt][i].mRstart == -1) continue;
            if(mSRs[svt][i].mRstart < (int32_t)hdr->target_len[mSRs[svt][i].mChr1]){
                mSRMapPos[mSRs[svt][i].mChr1].insert(std::make_pair(std::make_pair(mSRs[svt][i].mRstart, mSRs[svt][i].mID), mSRs[svt][i].mSVID));
            }
            if(mSRs[svt][i].mChr1 != mSRs[svt][i].mChr2 && mSRs[svt][i].mRstart < (int32_t)hdr->target_len[mSRs[svt][i].mChr2]){
                mSRMapPos[mSRs[svt][i].mChr2].insert(std::make_pair(std::make_pair(mSRs[svt][i].mRstart, mSRs[svt][i].mID), mSRs[svt][i].mSVID));
            }
        }
    }
    const uint16_t BAM_RDSKIP_MASK = (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FUNMAP);
    std::vector<std::multiset<std::string>> traSeqStore(svs.size()); // translocation SR read sequence
    std::vector<std::multiset<std::string>> triSeqStore(svs.size()); // translocation insertion sequence nearby bp
    std::vector<std::vector<uint8_t>> traQualStore(svs.size());      // translocation SR read mapping quality
    // Parse BAM
    for(auto& refIdx : mOpt->svRefID){
        if(mSRMapPos[refIdx].empty()) continue;
        // Load sequence
        int32_t seqlen = -1;
        char* chr1Seq = faidx_fetch_seq(fai, hdr->target_name[refIdx], 0, hdr->target_len[refIdx], &seqlen);
        std::vector<bool> hits(hdr->target_len[refIdx], false);
        // Collect all split-read pos 
        for(auto vsrIter = mSRMapPos[refIdx].begin(); vsrIter != mSRMapPos[refIdx].end(); ++vsrIter){
            hits[vsrIter->first.first] = true;
        }
        // Sequences and quality
        std::vector<std::multiset<std::string>> seqStore(svs.size());
        std::vector<std::multiset<std::string>> insStore(svs.size());
        std::vector<std::vector<uint8_t>> qualStore(svs.size());
        // Collect reads
        hts_itr_t* bamIter = sam_itr_queryi(idx, refIdx, 0, hdr->target_len[refIdx]);
        while(sam_itr_next(fp, bamIter, b) >= 0){
            if(b->core.flag & BAM_RDSKIP_MASK) continue;
            if(b->core.qual < mOpt->filterOpt->minMapQual || b->core.tid < 0) continue;
            if(!hits[b->core.pos]) continue;
            // Valid split-read
            size_t seed = svutil::hashString(bam_get_qname(b));
            auto vsrIter = mSRMapPos[refIdx].find(std::make_pair(b->core.pos, seed));
            if(vsrIter == mSRMapPos[refIdx].end()) continue;
            int32_t svid = vsrIter->second;
            int32_t svt = svs[svid].mSVT;
            // Get SR sequence
            std::string srseq = bamutil::getSeq(b);
            std::string siseq = "";
            int32_t bpInslen = 0;
            for(auto& rec : mSRs[svs[svid].mSVT]){
                if(rec.mID == seed){
                    bpInslen = rec.mInslen;
                    break;
                }
            }
            if(bpInslen > mOpt->filterOpt->minClipLen && svs[svid].mSVT != 4) svs[svid].getSCIns(b, srseq, siseq, bpInslen);
            // Adjust orientation
            bool bpPoint = false;
            if(svt >= 5){// translocation
                if(b->core.tid == svs[svid].mChr2) bpPoint = true;// bpPoint is true if b is on little chr
            }else{
                if(svt == 0){ // bpPoint is true if b is on 3' part of breakpoint
                    if(b->core.pos > svs[svid].mSVStart - mOpt->filterOpt->minClipLen) bpPoint = true;
                    else bpPoint = false;
                }else if(svs[svid].mSVT == 1){
                    if(b->core.pos > svs[svid].mSVEnd - mOpt->filterOpt->minClipLen) bpPoint = true;
                    else bpPoint = false;
                }
            }
            SRBamRecord::adjustOrientation(srseq, bpPoint, svt);
            if(!siseq.empty()) SRBamRecord::adjustOrientation(siseq, bpPoint, svt);
            // At most n split-reads used to to one SV event analysis
            if((int32_t)seqStore[svid].size() < mOpt->filterOpt->mMaxReadPerSV){
                if(svt >= 5){
                    traSeqStore[svid].insert(srseq);
                    traQualStore[svid].push_back(b->core.qual);
                    if(!siseq.empty()) triSeqStore[svid].insert(siseq);
                }else{
                    seqStore[svid].insert(srseq);
                    if(!siseq.empty()) insStore[svid].insert(siseq);
                    qualStore[svid].push_back(b->core.qual);
                }
            } 
        }
        // Process all SVs on this chromosome
        for(uint32_t svid = 0; svid < seqStore.size(); ++svid){
            if(svs[svid].mSVT >= 5) continue;
            if(svs[svid].mChr1 != refIdx) continue;
            // MSA
            bool bpRefined = false;
            if(seqStore[svid].size() > 1){
                AlignConfig alnCfg(5, -4, -10, -1, true, true);// both end gap free to keep each read ungapped as long as possible
                MSA* msa = new MSA(&seqStore[svid], mOpt->msaOpt->mMinCovForCS, mOpt->msaOpt->mMinBaseRateForCS, &alnCfg);
                msa->msa(svs[svid].mConsensus);
                delete msa;
                if(svs[svid].refineSRBp(mOpt, hdr, chr1Seq, chr1Seq)) bpRefined = true;
                if(!bpRefined){
                    svs[svid].mConsensus = "";
                    svs[svid].mSVRef = "";
                    svs[svid].mSRSupport = 0;
                    svs[svid].mSRAlignQuality = 0;
                    svs[svid].mSRMapQuality = 0;
                }else{// SR support and qualities
                    svs[svid].mSRSupport = seqStore[svid].size();
                    svs[svid].mSRMapQuality = statutil::median(qualStore[svid]);
                    if(insStore[svid].size() > 1){
                        MSA* imsa = new MSA(&insStore[svid], mOpt->msaOpt->mMinCovForCS, mOpt->msaOpt->mMinBaseRateForCS, &alnCfg);
                        imsa->msa(svs[svid].mBpInsSeq);
                        delete imsa;
                    }
                }
            }
        }
        if(chr1Seq) free(chr1Seq);
    }
    // Process translocations
    for(auto liteRefIdx = mOpt->svRefID.begin(); liteRefIdx != mOpt->svRefID.end(); ++liteRefIdx){
        char* liteChrSeq = NULL;
        auto largeRefIdx = liteRefIdx;
        ++largeRefIdx;
        for(; largeRefIdx !=mOpt->svRefID.end(); ++largeRefIdx){
            char* largeChrSeq = NULL;
            // Iterate SVs
            for(uint32_t svid = 0; svid < traSeqStore.size(); ++svid){
                if(svs[svid].mChr1 != (*largeRefIdx) || svs[svid].mChr2 != (*liteRefIdx)) continue;
                bool bpRefined = false;
                if(traSeqStore[svid].size() > 1){
                    if(!liteChrSeq){
                        int32_t liteChrSeqLen = -1;
                        liteChrSeq = faidx_fetch_seq(fai, hdr->target_name[*liteRefIdx], 0, hdr->target_len[*liteRefIdx], &liteChrSeqLen);
                    }
                    if(!largeChrSeq){
                        int32_t largeChrSeqLen = -1;
                        largeChrSeq = faidx_fetch_seq(fai, hdr->target_name[*largeRefIdx], 0, hdr->target_len[*largeRefIdx], &largeChrSeqLen);
                    }
                    AlignConfig alnCfg(5, -4, -10, -1, true, true);// both end gap free to keep each read ungapped as long as possible
                    MSA* msa = new MSA(&traSeqStore[svid], mOpt->msaOpt->mMinCovForCS, mOpt->msaOpt->mMinBaseRateForCS, &alnCfg);
                    msa->msa(svs[svid].mConsensus);
                    delete msa;
                    if(svs[svid].refineSRBp(mOpt, hdr, liteChrSeq, largeChrSeq)) bpRefined = true;
                    if(!bpRefined){
                        svs[svid].mConsensus = "";
                        svs[svid].mSVRef = "";
                        svs[svid].mSRSupport = 0;
                        svs[svid].mSRAlignQuality = 0;
                    }else{// SR support and qualities
                        svs[svid].mSRSupport = traSeqStore[svid].size();
                        svs[svid].mSRMapQuality = util::median(traQualStore[svid]);
                        if(triSeqStore.size() > 0){
                            MSA* imsa = new MSA(&triSeqStore[svid], mOpt->msaOpt->mMinCovForCS, mOpt->msaOpt->mMinBaseRateForCS, &alnCfg);
                            imsa->msa(svs[svid].mBpInsSeq);
                            delete imsa;
                        }
                    }
                }
            }
            if(largeChrSeq) free(largeChrSeq);
        }
        if(liteChrSeq) free(liteChrSeq);
    }
    // Add ChrName
    for(uint32_t i = 0; i < svs.size(); ++i){
        svs[i].mNameChr1 = hdr->target_name[svs[i].mChr1];
        svs[i].mNameChr2 = hdr->target_name[svs[i].mChr2];
    }
    // Clean-up
    sam_close(fp);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    bam_destroy1(b);
}
