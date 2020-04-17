#include "srbamrecord.h"
#include "breakpoint.h"

void SRBamRecordSet::classifyJunctions(JunctionMap* jctMap){
    util::loginfo("Beg classifing SRs into various SV candidates");
    int svtIdx = 0;
    int32_t rst = -1;
    bool srchr1 = false;
    uint32_t j = 0;
    int32_t inslen = 0;
    for(auto iter = jctMap->mJunctionReads.begin(); iter != jctMap->mJunctionReads.end(); ++iter){
        for(uint32_t i = 0; i < iter->second.size(); i += 2){
            j = i + 1;
            if(iter->second[i].mSCLen && iter->second[j].mSCLen && 
               (iter->second[j].mSCLen < iter->second[i].mSeqmatch - mOpt->filterOpt->mMaxReadSep ||
               iter->second[i].mSCLen < iter->second[j].mSeqmatch - mOpt->filterOpt->mMaxReadSep)){
                continue;
            }
            if(iter->second[i].mSCLen > iter->second[j].mSeqmatch){
                inslen = iter->second[i].mSCLen - iter->second[j].mSeqmatch;
            }else{
                inslen = 0;
            }
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
                    rst = iter->second[largerChrIdx].mRstart;
                    srchr1 = true;
                    if(rst < 0){
                        rst = iter->second[littleChrIdx].mRstart;
                        srchr1 = false;
                    }
                    mSRs[svtIdx].push_back(SRBamRecord(iter->second[largerChrIdx].mRefidx,
                                                    iter->second[largerChrIdx].mRefpos,
                                                    iter->second[littleChrIdx].mRefidx,
                                                    iter->second[littleChrIdx].mRefpos,
                                                    rst,
                                                    srchr1,
                                                    inslen,
                                                    iter->first,
                                                    iter->second[littleChrIdx].mRead1));
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
                    if(iter->second[leftPart].mSCleft) std::swap(leftPart, rightPart);
                    svtIdx = 4; 
                }else if(iter->second[j].mForward == iter->second[i].mForward && // same direction
                         iter->second[j].mSCleft != iter->second[i].mSCleft && // opposing soft-clips
                         std::abs(iter->second[j].mRefpos - iter->second[i].mRefpos) >= mOpt->filterOpt->mMinRefSep){// breakpoint faraway
                    if(iter->second[leftPart].mSCleft) svtIdx = 3; // left part leading soft-clip, duplication
                    else svtIdx = 2; // left part tailing soft-clip, deletion
                }else if(iter->second[j].mForward != iter->second[i].mForward && // opposing direction
                         iter->second[j].mSCleft == iter->second[i].mSCleft && // same soft-clips
                         std::abs(iter->second[j].mRefpos - iter->second[i].mRefpos) >= mOpt->filterOpt->mMinRefSep){// breakpoint farway
                    if(iter->second[rightPart].mSCleft) svtIdx = 1; // 3to3 right spanning inversion breakpoint
                    else svtIdx = 0; // 5to5 left spanning inversion breakpoint
                }
                if(svtIdx != -1 && mOpt->SVTSet.find(svtIdx) != mOpt->SVTSet.end()){
                    rst = iter->second[leftPart].mRstart;
                    srchr1 = true;
                    if(rst < 0){
                        rst = iter->second[rightPart].mRstart;
                        srchr1 = false;
                    }
                    mSRs[svtIdx].push_back(SRBamRecord(iter->second[leftPart].mRefidx,
                                                       iter->second[leftPart].mRefpos,
                                                       iter->second[rightPart].mRefidx,
                                                       iter->second[rightPart].mRefpos,
                                                       rst,
                                                       srchr1,
                                                       inslen,
                                                       iter->first,
                                                       iter->second[leftPart].mRead1));
                }
            }
        }
    }
    util::loginfo("End classifing SRs into various SV candidates");
}

void SRBamRecordSet::cluster(std::vector<SRBamRecord>& srs, SVSet& svs, int32_t svt){
    int32_t origSize = svs.size();
    util::loginfo("Beg clustering SRs for SV type " + std::to_string(svt) + ", all " + std::to_string(srs.size()) + " SRs ");
#ifdef DEBUG
    if(mOpt->debug & DEBUG_FCALL){
        std::cout << "DEBUG_BEG_CLUSTER_SRS_FOR_SV_TYPE:" << svt << std::endl;
    }
#endif
    std::set<int32_t> clique; // component cluster
    int32_t totsrs = srs.size();
    int32_t i = 0, j = 0;
    while(i < totsrs){
        clique.clear();
        clique.insert(i);
        j = i + 1;
        while(j < totsrs){
            if(srs[j].mChr1 == srs[i].mChr1 && 
               srs[j].mChr2 == srs[i].mChr2 &&
               srs[j].mPos1 - srs[i].mPos1 < mOpt->filterOpt->mMaxReadSep &&
               std::abs(srs[j].mPos2 - srs[i].mPos2) < mOpt->filterOpt->mMaxReadSep){
                clique.insert(j);
                ++j;
            }else{
                break;
            }
        }
        searchCliques(clique, srs, svs, svt);
        i = j;
    }
    util::loginfo("End clustering SRs for SV type " + std::to_string(svt) + ", got " + std::to_string(svs.size() - origSize) + " SV candidates.");
#ifdef DEBUG
    if(mOpt->debug & DEBUG_FCALL){
        std::cout << "End clustering SRs for SV type " << svt << std::endl;
    }
#endif
}

void SRBamRecordSet::searchCliques(std::set<int32_t>& clique, std::vector<SRBamRecord>& srs, SVSet& svs, int32_t svt){
    auto iter = clique.begin();
    int32_t srid = *iter;
    int32_t ciposlow = srs[srid].mPos1;
    int32_t ciposhigh = srs[srid].mPos1;
    int32_t ciendlow = srs[srid].mPos2;
    int32_t ciendhigh = srs[srid].mPos2;
    uint64_t pos1 = srs[srid].mPos1;
    uint64_t pos2 = srs[srid].mPos2;
    int32_t inslen = 0, inscnt = 0;
    if(srs[srid].mInslen >= mOpt->filterOpt->mMinBpInsLen){
        inslen += srs[srid].mInslen;
        inscnt = 1;
    }
    int32_t chr1 = srs[srid].mChr1;
    int32_t chr2 = srs[srid].mChr2;
    ++iter;
    for(; iter != clique.end(); ++iter){
        srid = *iter;
        ciposlow = std::min(srs[srid].mPos1, ciposlow);
        ciposhigh = std::max(srs[srid].mPos1, ciposhigh);
        ciendlow = std::min(srs[srid].mPos2, ciendlow);
        ciendhigh = std::max(srs[srid].mPos2, ciendhigh);
        pos1 += srs[srid].mPos1;
        pos2 += srs[srid].mPos2;
        if(srs[srid].mInslen >= mOpt->filterOpt->mMinBpInsLen){
            inslen += srs[srid].mInslen;
            ++inscnt;
        }
    }
    if(clique.size() == 1 && 
       inscnt == 1 && 
       inslen > mOpt->filterOpt->mMaxSingSrSeedIns){ // insertion length too long singleton sr seed, drop
        return;
    }

    if(clique.size() >= mOpt->filterOpt->mMinSeedSR){
        int32_t svStart = pos1/clique.size();
        int32_t svEnd = pos2/clique.size();
        bool hotfusion = false;
        if(mOpt->pairOlpRegs &&
           mOpt->pairOlpRegs->overlap(sam_hdr_tid2name(mOpt->bamheader, chr1), svStart, svStart + 1) &&
           mOpt->pairOlpRegs->overlap(sam_hdr_tid2name(mOpt->bamheader, chr2), svEnd, svEnd + 1)){
            hotfusion = true;
        }
        if((hotfusion && (int32_t)clique.size() >= mOpt->fuseOpt->mWhiteFilter.mMinSRSeed) ||
           ((!hotfusion) && (int32_t)clique.size() >= mOpt->fuseOpt->mUsualFilter.mMinSRSeed)){
            int32_t svISize = 0;
            int32_t svIScnt = 0;
            if(inscnt > .3 * clique.size()){
                svISize = inslen / inscnt; 
                svIScnt = inscnt;
            }
            int32_t svid = svs.size();
            SVRecord* svr = new SVRecord();
            svr->mChr1 = chr1;
            svr->mSVStart = svStart;
            svr->mChr2 = chr2;
            svr->mSVEnd = svEnd;
            svr->mCiPosLow = ciposlow - svStart;
            svr->mCiPosHigh = ciposhigh - svStart;
            svr->mCiEndLow = ciendlow - svEnd;
            svr->mCiEndHigh = ciendhigh - svEnd;
            svr->mAlnInsLen = svISize;
            svr->mInsSeedCnt = svIScnt;
            svr->mID = svid;
            svr->mSVT = svt;
            if(clique.size() == 1) svr->mFromOneSR = true;
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

void SRBamRecordSet::assembleOneContig(SVSet& svs, int32_t refIdx){
    if(mSRMapPos[refIdx].empty()){
        if(mOpt->svRefID.find(refIdx) != mOpt->svRefID.end()){
            faidx_t* fai = fai_load(mOpt->alnref.c_str());
            int32_t seqlen = -1;
            char* chr1Seq = faidx_fetch_seq(fai, mOpt->bamheader->target_name[refIdx], 0, mOpt->bamheader->target_len[refIdx], &seqlen);
            for(uint32_t svid = 0; svid < svs.size(); ++svid){
                if(svs[svid]->mChr1 != refIdx && svs[svid]->mChr2 != refIdx) continue;
                if(svs[svid]->mSVT >= 5 && svs[svid]->mChr2 == refIdx){// Fetch lite chr seq of translocation
                    BreakPoint bp = BreakPoint(svs[svid], mOpt->bamheader, 10 * mOpt->libInfo->mReadLen);
                    svs[svid]->mTraChr2Seq = bp.getLittleRef(chr1Seq);
                    continue;
                }
                if(svs[svid]->mSVT >= 5 && svs[svid]->mChr1 == refIdx){// Fetch large chr seq of translocation
                    BreakPoint bp = BreakPoint(svs[svid], mOpt->bamheader, 10 * mOpt->libInfo->mReadLen);
                    svs[svid]->mTraChr1Seq = bp.getLargerRef(chr1Seq);
                    continue;
                }
            }
            if(chr1Seq){
                free(chr1Seq);
                chr1Seq = NULL;
            }
            fai_destroy(fai);
            fai = NULL;
            return;
        }
    }
    // Open file handles
    samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
    hts_set_fai_filename(fp, mOpt->alnref.c_str());
    hts_idx_t* idx = sam_index_load(fp, mOpt->bamfile.c_str());
    faidx_t* fai = fai_load(mOpt->alnref.c_str());
    bam1_t* b = bam_init1();
    util::loginfo("Beg assembling SRs on contig: " + std::string(mOpt->bamheader->target_name[refIdx]), mOpt->logMtx);
    const uint16_t BAM_RDSKIP_MASK = (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FUNMAP | BAM_FSUPPLEMENTARY);
    // Load sequence
    int32_t seqlen = -1;
    char* chr1Seq = faidx_fetch_seq(fai, mOpt->bamheader->target_name[refIdx], 0, mOpt->bamheader->target_len[refIdx], &seqlen);
    std::vector<bool> hits(mOpt->bamheader->target_len[refIdx], false);
    // Collect all split-read pos
    for(auto vsrIter = mSRMapPos[refIdx].begin(); vsrIter != mSRMapPos[refIdx].end(); ++vsrIter){
        hits[vsrIter->first.first] = true;
    }
    // Sequences and quality
    std::vector<std::multiset<std::string>> seqStore(svs.size());
    std::vector<std::multiset<std::string>> insStore(svs.size());
    std::vector<std::vector<uint8_t>> qualStore(svs.size());
    // Collect reads
    hts_itr_t* bamIter = sam_itr_queryi(idx, refIdx, 0, mOpt->bamheader->target_len[refIdx]);
    while(sam_itr_next(fp, bamIter, b) >= 0){
        if(b->core.flag & BAM_RDSKIP_MASK) continue;
        if(b->core.qual < mOpt->filterOpt->minMapQual || b->core.tid < 0) continue;
        if(!hits[b->core.pos]) continue;
        // Valid split-read
        size_t seed = svutil::hashString(bam_get_qname(b));
        auto vsrIter = mSRMapPos[refIdx].find(std::make_pair(b->core.pos, seed));
        if(vsrIter == mSRMapPos[refIdx].end()) continue;
        int32_t svid = vsrIter->second;
        int32_t svt = svs[svid]->mSVT;
        int32_t bpInslen = 0;
        bool r12mismatch = true;
        bool oread1 = false;
        for(auto& rec : mSRs[svs[svid]->mSVT]){
            oread1 = (b->core.flag & BAM_FREAD1);
            if(rec.mID == seed && rec.mRead1 == oread1){
                bpInslen = rec.mInslen;
                r12mismatch = false;
                break;
            }
        }
        if(r12mismatch) continue;
        std::string srseq = ""; // read sequence excluding insertiong after clip pos
        std::string siseq = ""; // inserted sequence after clip pos
        svs[svid]->getSCIns(b, srseq, siseq, bpInslen, mOpt->filterOpt->mMinBpInsLen);
        // Adjust orientation
        bool bpPoint = false;
        if(svt >= 5 && b->core.tid == svs[svid]->mChr2) bpPoint = true;
        else{
            if(svt == 0){
                if(svs[svid]->mSVStart - b->core.pos < mOpt->filterOpt->minClipLen) bpPoint = true;
                else bpPoint = false;
            }else if(svt == 1){
                if(svs[svid]->mSVEnd - b->core.pos < mOpt->filterOpt->minClipLen) bpPoint = true;
                else bpPoint = false;
            }
        }
        SRBamRecord::adjustOrientation(srseq, bpPoint, svt);
        if(!siseq.empty()) SRBamRecord::adjustOrientation(siseq, bpPoint, svt);
        // At most n split-reads used to to one SV event analysis
        if(svt >= 5){
            mLock.lock();
            if((int32_t)mTraSeqStore[svid].size() < mOpt->filterOpt->mMaxReadPerSV){
                if(!srseq.empty()){
                    mTraSeqStore[svid].insert(srseq);
                    mTraQualStore[svid].push_back(b->core.qual);
                }
                if(!siseq.empty()){
                    mTriSeqStore[svid].insert(siseq);
                }
            }
            mLock.unlock();
        }else{
            if((int32_t)seqStore[svid].size() < mOpt->filterOpt->mMaxReadPerSV){
                if(!srseq.empty()){
                    seqStore[svid].insert(srseq);
                    qualStore[svid].push_back(b->core.qual);
                }
                if(!siseq.empty()){
                    insStore[svid].insert(siseq);
                }
            }
        } 
    }
    // Process all SVs on this chromosome
    for(uint32_t svid = 0; svid < seqStore.size(); ++svid){
        if(svs[svid]->mChr1 != refIdx && svs[svid]->mChr2 != refIdx) continue;
        if(svs[svid]->mSVT >= 5 && svs[svid]->mChr2 == refIdx){// Fetch lite chr seq of translocation
            BreakPoint bp = BreakPoint(svs[svid], mOpt->bamheader, 10 * mOpt->libInfo->mReadLen);
            svs[svid]->mTraChr2Seq = bp.getLittleRef(chr1Seq);
            continue;
        }
        if(svs[svid]->mSVT >= 5 && svs[svid]->mChr1 == refIdx){// Fetch large chr seq of translocation
            BreakPoint bp = BreakPoint(svs[svid], mOpt->bamheader, 10 * mOpt->libInfo->mReadLen);
            svs[svid]->mTraChr1Seq = bp.getLargerRef(chr1Seq);
            continue;
        }
        // MSA
        bool bpRefined = false;
        if(seqStore[svid].size()){
            if(svs[svid]->mFromOneSR){
                svs[svid]->mConsensus = *(seqStore[svid].begin());
                if(insStore[svid].size()) svs[svid]->mBpInsSeq = *(insStore[svid].begin());
            }else{
                AlignConfig alnCfg(5, -4, -10, -1, true, true);// both end gap free to keep each read ungapped as long as possible
                MSA* msa = new MSA(&seqStore[svid], mOpt->msaOpt->mMinCovForCS, mOpt->msaOpt->mMinBaseRateForCS, &alnCfg);
                if(seqStore[svid].size() < 3) msa->mMinCovForCS = seqStore[svid].size();
                msa->msa(svs[svid]->mConsensus);
                delete msa;
                if(insStore[svid].size()){
                    if(insStore[svid].size() == 1){
                        svs[svid]->mBpInsSeq = *(insStore[svid].begin());
                    }else{
                        MSA* imsa = new MSA(&insStore[svid], mOpt->msaOpt->mMinCovForCS, mOpt->msaOpt->mMinBaseRateForCS, &alnCfg);
                        if(insStore[svid].size() < 3) imsa->mMinCovForCS = insStore[svid].size();
                        imsa->msa(svs[svid]->mBpInsSeq);
                        delete imsa;
                    }
                }
            }
            if(svs[svid]->refineSRBp(mOpt, mOpt->bamheader, chr1Seq, chr1Seq)) bpRefined = true;
            if(!bpRefined){
#ifdef DEBUG
                if(mOpt->debug & DEBUG_FCALL){
                    std::cout << "DEBUG_BPREFINED_FAILED_SV:\n" << svs[svid] << std::endl;
                }
#endif
                svs[svid]->mConsensus = "";
                svs[svid]->mSVRef = "";
                svs[svid]->mSRSupport = 0;
                svs[svid]->mSRAlignQuality = 0;
                svs[svid]->mSRMapQuality = 0;
            }else{// SR support and qualities
                svs[svid]->mSVRef = svs[svid]->mSVRef.substr(svs[svid]->mGapCoord[2] - 1, 1);
                svs[svid]->mSRSupport = seqStore[svid].size();
                svs[svid]->mSRMapQuality = statutil::median(qualStore[svid]);
                svs[svid]->mBpInsFreq = (float)(insStore[svid].size())/(float)(seqStore[svid].size());
            }
            svs[svid]->mTraChr1Seq.clear(); svs[svid]->mTraChr1Seq.shrink_to_fit();
            svs[svid]->mTraChr2Seq.clear(); svs[svid]->mTraChr2Seq.shrink_to_fit();
        }
    }
    if(chr1Seq) free(chr1Seq);
    util::loginfo("End assembling SRs on contig: " + std::string(mOpt->bamheader->target_name[refIdx]), mOpt->logMtx);
    sam_close(fp);
    hts_idx_destroy(idx);
    fai_destroy(fai);
    bam_destroy1(b);
}

void SRBamRecordSet::assembleSplitReads(SVSet& svs){
    // Construct bool filter of SR mapping positions
    for(uint32_t svt = 0; svt < mSRs.size(); ++svt){
        for(uint32_t i = 0; i < mSRs[svt].size(); ++i){
            if(mSRs[svt][i].mSVID == -1) continue;
            if(mSRs[svt][i].mSRChr1){
                mSRMapPos[mSRs[svt][i].mChr1].insert(std::make_pair(std::make_pair(mSRs[svt][i].mRstart, mSRs[svt][i].mID), mSRs[svt][i].mSVID));
            }else{
                mSRMapPos[mSRs[svt][i].mChr2].insert(std::make_pair(std::make_pair(mSRs[svt][i].mRstart, mSRs[svt][i].mID), mSRs[svt][i].mSVID));
            }
        }
    }
    // reserve space for translocation information
    mTraSeqStore.resize(svs.size());
    mTraQualStore.resize(svs.size());
    mTriSeqStore.resize(svs.size());
    // parallel assemble on each contig
    std::vector<std::future<void>> asret(mOpt->svRefID.size());
    int i = 0;
    for(auto& refIdx: mOpt->svRefID) asret[i++] = mOpt->pool->enqueue(&SRBamRecordSet::assembleOneContig, this, std::ref(svs), refIdx);
    for(auto& e: asret) e.get();
    // assemble translocation srs
    util::loginfo("Beg assembling SRs across chromosomes", mOpt->logMtx);
    std::vector<int32_t> nonOnesIds;
    std::vector<int32_t> yesOnesIds;
    for(uint32_t cri = 0; cri < svs.size(); ++cri){
        if(svs[cri]->mSVT >= 5){
            if(svs[cri]->mFromOneSR) yesOnesIds.push_back(cri);
            else nonOnesIds.push_back(cri);
        }
    }
    std::vector<std::pair<int32_t, int32_t>> nspidx, yspidx;
    int32_t nps = util::divideVecIdx(nonOnesIds.size(), mOpt->nthread, nspidx);
    int32_t yps = util::divideVecIdx(yesOnesIds.size(), mOpt->nthread, yspidx);
    std::vector<int32_t> crsvids;
    int32_t cps = 0;
    for(cps = 0; cps < std::min(nps, yps); ++cps){
        for(int32_t insp = nspidx[cps].first; insp < nspidx[cps].second; ++insp){
            crsvids.push_back(nonOnesIds[insp]);
        }
        for(int32_t iysp = yspidx[cps].first; iysp < yspidx[cps].second; ++iysp){
            crsvids.push_back(yesOnesIds[iysp]);
        }
    }
    if(cps < nps){
        for(; cps < nps; ++cps){
            for(int32_t insp = nspidx[cps].first; insp < nspidx[cps].second; ++insp){
                crsvids.push_back(nonOnesIds[insp]);
            }
        }
    }
    if(cps < yps){
        for(; cps < yps; ++cps){
            for(int32_t iysp = yspidx[cps].first; iysp < yspidx[cps].second; ++iysp){
                crsvids.push_back(yesOnesIds[iysp]);
            }
        }
    }
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    util::divideVecIdx(crsvids.size(), mOpt->nthread, vpidx);
    AlignConfig* alnCfg = new AlignConfig(5, -4, -10, -1, true, true);
    std::vector<std::future<void>> casret(vpidx.size());
    for(uint32_t pidx = 0; pidx < vpidx.size(); ++pidx){
        casret[pidx] = mOpt->pool->enqueue(&SRBamRecordSet::assembleCrossChr, this, std::ref(svs), alnCfg, std::ref(crsvids), vpidx[pidx].first, vpidx[pidx].second);
    }
    for(auto& e: casret) e.get();
    // Add ChrName
    for(uint32_t i = 0; i < svs.size(); ++i){
        svs[i]->mNameChr1 = mOpt->bamheader->target_name[svs[i]->mChr1];
        svs[i]->mNameChr2 = mOpt->bamheader->target_name[svs[i]->mChr2];
    }
    // Clean-up
    delete alnCfg;
    util::loginfo("End assembling SRs across chromosomes", mOpt->logMtx);
}

void SRBamRecordSet::assembleCrossChr(SVSet& svs, AlignConfig* alnCfg, const std::vector<int32_t>& crsidx, int32_t begIdx, int32_t endIdx){
    util::loginfo("Beg processing vpidx: " + std::to_string(begIdx) + "->" + std::to_string(endIdx), mOpt->logMtx);
    for(int32_t cidx = begIdx; cidx < endIdx; ++cidx){
        int32_t svid = crsidx[cidx];
        bool bpRefined = false;
        if(mTraSeqStore[svid].size()){
            if(svs[svid]->mFromOneSR){
                svs[svid]->mConsensus = *(mTraSeqStore[svid].begin());
                if(mTriSeqStore[svid].size()) svs[svid]->mBpInsSeq = *(mTriSeqStore[svid].begin());
            }else{
                MSA* msa = new MSA(&mTraSeqStore[svid], mOpt->msaOpt->mMinCovForCS, mOpt->msaOpt->mMinBaseRateForCS, alnCfg);
                if(mTraSeqStore[svid].size() < 3) msa->mMinCovForCS = mTraSeqStore[svid].size();
                msa->msa(svs[svid]->mConsensus);
                delete msa;
                if(mTriSeqStore[svid].size()){
                    if(mTriSeqStore[svid].size() == 1){
                        svs[svid]->mBpInsSeq = *(mTriSeqStore[svid].begin());
                    }else{
                        MSA* imsa = new MSA(&mTriSeqStore[svid], mOpt->msaOpt->mMinCovForCS, mOpt->msaOpt->mMinBaseRateForCS, alnCfg);
                        if(mTriSeqStore[svid].size() < 3) imsa->mMinCovForCS = mTriSeqStore[svid].size();
                        imsa->msa(svs[svid]->mBpInsSeq);
                        delete imsa;
                    }
                }
            }
            if(svs[svid]->refineSRBp(mOpt, mOpt->bamheader, NULL, NULL)) bpRefined = true;
            if(!bpRefined){
#ifdef DEBUG
                if(mOpt->debug & DEBUG_FCALL){
                    std::cout << "DEBUG_BPREFINED_FAILED_SV:\n" << svs[svid] << std::endl;
                }
#endif
                svs[svid]->mConsensus = "";
                svs[svid]->mSVRef = "";
                svs[svid]->mSRSupport = 0;
                svs[svid]->mSRAlignQuality = 0;
            }else{// SR support and qualities
                svs[svid]->mSVRef = svs[svid]->mSVRef.substr(svs[svid]->mGapCoord[2] - 1, 1);
                svs[svid]->mSRSupport = mTraSeqStore[svid].size();
                svs[svid]->mSRMapQuality = util::median(mTraQualStore[svid]);
                svs[svid]->mBpInsFreq = (float)(mTriSeqStore[svid].size())/(float)(mTraSeqStore[svid].size());
            }
            svs[svid]->mTraChr1Seq.clear(); svs[svid]->mTraChr1Seq.shrink_to_fit();
            svs[svid]->mTraChr2Seq.clear(); svs[svid]->mTraChr2Seq.shrink_to_fit();
        }
    }
    util::loginfo("End processing vpidx: " + std::to_string(begIdx) + "->" + std::to_string(endIdx), mOpt->logMtx);
}
