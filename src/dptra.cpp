#include "stats.h"

TraDPStat* TraDPStat::merge(const std::vector<TraDPStat*>& sts, int32_t n){
    TraDPStat* ret = new TraDPStat(n);
    for(int32_t j = 0; j < n; ++j){
        for(uint32_t i = 0; i < sts.size(); ++i){
            ret->mSpnCnts[j].mAlth1 += sts[i]->mSpnCnts[j].mAlth1;
            ret->mSpnCnts[j].mAlth2 += sts[i]->mSpnCnts[j].mAlth2;
            ret->mSpnCnts[j].mAltQual.insert(ret->mSpnCnts[j].mAltQual.end(), sts[i]->mSpnCnts[j].mAltQual.begin(), sts[i]->mSpnCnts[j].mAltQual.end());
        }
    }
    return ret;
}

void TraDPStat::dptra(const ContigSpanPoints& spPts){
    samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    util::loginfo("Start processing alt DPs on: " + std::string(h->target_name[mBigChr]) + " and " + std::string(h->target_name[mLiteChr]), mOpt->logMtx);
    hts_idx_t* idx = sam_index_load(fp, mOpt->bamfile.c_str());
    // Flag spanning breakpoints
    std::set<int32_t> spnBp;
    for(uint32_t i = 0; i < spPts[mBigChr].size(); ++i) spnBp.insert(spPts[mBigChr][i].mBpPos);
    // Count reads on mLiteChr
    hts_itr_t* itr = sam_itr_queryi(idx, mLiteChr, 0, h->target_len[mLiteChr]);
    bam1_t* b = bam_init1();
    const uint16_t COV_STAT_SKIP_MASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FMUNMAP);
    std::unordered_map<size_t, uint8_t> qualities;
    while(sam_itr_next(fp, itr, b) >= 0){
        if(!(b->core.flag & BAM_FPAIRED)) continue;
        if(b->core.flag & COV_STAT_SKIP_MASK) continue;
        if(b->core.qual < mOpt->filterOpt->mMinGenoQual) continue;
        if(b->core.mtid != mBigChr) continue;
        size_t hv = svutil::hashPairCurr(b);
        qualities[hv] = b->core.qual;
    }
    hts_itr_destroy(itr);
    // Count reads on mLiteChr
    itr = sam_itr_queryi(idx, mBigChr, 0, h->target_len[mBigChr]);
    while(sam_itr_next(fp, itr, b) >= 0){
        if(!(b->core.flag & BAM_FPAIRED)) continue;
        if(b->core.flag & COV_STAT_SKIP_MASK) continue;
        if(b->core.qual < mOpt->filterOpt->mMinGenoQual) continue;
        if(b->core.mtid != mLiteChr) continue;
        size_t hv = svutil::hashPairMate(b);
        auto itq = qualities.find(hv);
        if(itq == qualities.end()) continue;
        uint8_t pairQual = std::min(itq->second, b->core.qual);
        itq->second = 0;
        // Pair quality
        if(pairQual < mOpt->filterOpt->mMinGenoQual) continue; // Low quality pair
        // Get SV type
        int32_t svt =  DPBamRecord::getSVType(b, mOpt);
        if(svt == -1) continue;
        // Spanning a breakpoint?
        bool spanvalid = false;
        int32_t pbegin = b->core.pos;
        int32_t pend = std::min(b->core.pos + mOpt->libInfo->mMaxNormalISize, (int32_t)h->target_len[mBigChr]);
        if(b->core.flag & BAM_FREVERSE){
            pbegin = std::max(0, b->core.pos + b->core.l_qseq - mOpt->libInfo->mMaxNormalISize);
            pend = std::min(b->core.pos + b->core.l_qseq, (int32_t)h->target_len[mBigChr]);
        }
        for(int32_t i = pbegin; i < pend; ++i){
            if(spnBp.find(i) != spnBp.end()){
                spanvalid = true;
                break;
            }
        }
        if(spanvalid){
            // Fetch all relevant SVs
            auto itspan = std::lower_bound(spPts[mBigChr].begin(), spPts[mBigChr].end(), SpanPoint(pbegin));
            for(; itspan != spPts[mBigChr].end() && pend >= itspan->mBpPos; ++itspan){
                if(svt == itspan->mSVT){
                    uint8_t* hpptr = bam_aux_get(b, "HP");
                    mSpnCnts[itspan->mID].mAltQual.push_back(pairQual);
                    if(hpptr){
                        mOpt->libInfo->mIsHaploTagged = true;
                        int hap = bam_aux2i(hpptr);
                        if(hap == 1) ++mSpnCnts[itspan->mID].mAlth1;
                        else ++mSpnCnts[itspan->mID].mAlth2;
                    }
                    mOpt->outMtx.lock();
                    bam_aux_update_int(b, "SVID", itspan->mID);
                    assert(sam_write1(mOpt->fbamout, h, b) >= 0);
                    mOpt->outMtx.unlock();
                }
            }
        }
    }
    util::loginfo("Finish processing alt DPs on: " + std::string(h->target_name[mBigChr]) + " and " + std::string(h->target_name[mLiteChr]), mOpt->logMtx);
    // Clean-up
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
}
