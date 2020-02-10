#include "svrecord.h"
#include "breakpoint.h"

template<typename TBreakPoint>
bool SVRecord::coordTransform(TBreakPoint& bp, AlignDescriptor& ad, int32_t& finalGapStart, int32_t& finalGapEnd){
    if(mSVT == 0 || mSVT == 5){// 5to5
        int32_t annealed = bp.mSVStartEnd - bp.mSVStartBeg;
        if(ad.mRefStart >= annealed || (ad.mRefEnd < annealed)) return false;
        finalGapStart = bp.mSVStartBeg + ad.mRefStart;
        finalGapEnd = bp.mSVEndBeg + mSVRef.size() - ad.mRefEnd;
        return true;
    }
    if(mSVT == 1 || mSVT == 6){// 3to3
        int32_t annealed = bp.mSVStartEnd - bp.mSVStartBeg;
        if(ad.mRefStart >= annealed || (ad.mRefEnd < annealed)) return false;
        finalGapStart = bp.mSVStartBeg + (annealed - ad.mRefStart);
        finalGapEnd = bp.mSVEndBeg + (ad.mRefEnd - annealed);
        return true;
    }
    if(mSVT == 2 || mSVT == 7){// 5to3
        int32_t annealed = bp.mSVStartEnd - bp.mSVStartBeg;
        if(ad.mRefStart >= annealed || (ad.mRefEnd < annealed)) return false;
        finalGapStart = bp.mSVStartBeg + ad.mRefStart;
        finalGapEnd = bp.mSVEndBeg + (ad.mRefEnd - annealed);
        return true;
    }
    if(mSVT == 3 || mSVT == 8){// 3to5
        int32_t annealed = bp.mSVEndEnd - bp.mSVEndBeg;
        if(ad.mRefStart >= annealed || (ad.mRefEnd < annealed)) return false;
        finalGapStart = bp.mSVStartBeg + (ad.mRefEnd - annealed);
        finalGapEnd = bp.mSVEndBeg + ad.mRefStart;
        return true;
    }
    if(mSVT == 4){// insertion
        finalGapStart = bp.mSVStartBeg + ad.mRefStart;
        finalGapEnd = bp.mSVStartBeg + ad.mRefStart;
        return true;
    }
    return false;
}

void SVRecord::findHomology(AlignDescriptor& ad){
    if(mSVT == 4){
        ad.mHomRight = Aligner::longestHomology(mConsensus.substr(ad.mCSStart), mSVRef.substr(ad.mRefEnd - 1), 1);
        std::string preC = mConsensus.substr(0, ad.mCSEnd - 1);
        std::string preR = mSVRef.substr(0, ad.mRefStart);
        util::reverse(preC);
        util::reverse(preR);
        ad.mHomLeft = Aligner::longestHomology(preC, preR, 1);
    }else{
        ad.mHomRight = Aligner::longestHomology(mConsensus.substr(ad.mCSEnd - 1), mSVRef.substr(ad.mRefStart), 1);
        std::string preC = mConsensus.substr(0, ad.mCSStart);
        std::string preR = mSVRef.substr(0, ad.mRefEnd - 1);
        util::reverse(preC);
        util::reverse(preR);
        ad.mHomLeft = Aligner::longestHomology(preC, preR, 1);
    }
}

bool SVRecord::findSplit(Matrix2D<char>* alnResult, AlignDescriptor& ad){
    // Initialization
    int32_t gs = 0; // gap start
    int32_t ge = 0; // gap end
    // Find longest internal gap
    int32_t refIndex = 0;
    int32_t varIndex = 0;
    int32_t gapStartRefIndex = 0;
    int32_t gapStartVarIndex = 0;
    int32_t a1 = 0;
    bool inGap = false;
    for(int32_t j = 0; j < alnResult->ncol(); ++j){
        if(alnResult->get(0, j) != '-') ++varIndex;
        if(alnResult->get(1, j) != '-') ++refIndex;
        // Internal gap status check
        if((alnResult->get(0, j) == '-' || alnResult->get(1, j) == '-') && varIndex && refIndex){
            if(!inGap){
                gapStartVarIndex = (alnResult->get(0, j) != '-') ? varIndex - 1 : varIndex;
                gapStartRefIndex = (alnResult->get(1, j) != '-') ? refIndex - 1 : refIndex;
                a1 = j;
                inGap = true;
            }
        }else{
            if(inGap && checkSVGap(refIndex - gapStartRefIndex, ad.mRefEnd - ad.mRefStart, varIndex - gapStartVarIndex, ad.mCSEnd - ad.mCSStart)){
                ad.mRefStart = gapStartRefIndex;
                ad.mRefEnd = refIndex;
                ad.mCSStart = gapStartVarIndex;
                ad.mCSEnd = varIndex;
                gs = a1;
                ge = j - 1;
            }
            inGap = false;
        }
    }
    if(ad.mRefEnd < ad.mRefStart) return false;
    // check whether this is a valid split-read alignment
    if(!ad.validGapAlignment(mSVT)) return false;
    // check percent identity
    ad.mPercID = alnResult->identityPercent(gs, ge);
    if(!ad.validFlankQual()) return false;
    // find homology
    findHomology(ad);
    // Valid split-read alignment
    return true;
} 

    
bool SVRecord::consensusRefAlign(Matrix2D<char>* alnResult){
    AlignConfig* alnCfg = new AlignConfig(5, -4, -4, -4, true, false);
    Aligner* aligner = new Aligner();
    aligner->mAlignConfig = alnCfg;
    if(mSVT == 4){
        bool alnRet = aligner->splitAligner(mSVRef, mConsensus, alnResult);
        for(int j = 0; j < alnResult->ncol(); ++j){
            char tmp = alnResult->get(0, j);
            alnResult->set(0, j) = alnResult->get(1, j);
            alnResult->set(1, j) = tmp;
        }
        delete alnCfg;
        delete aligner;
        return alnRet;
    }else{
        alnCfg->mHorizontalEndGapFree = true;
        alnCfg->mVerticalEndGapFree = false;
        bool alnRet = aligner->splitAligner(mConsensus, mSVRef, alnResult);
        delete alnCfg;
        delete aligner;
        return alnRet;
    }
}

bool SVRecord::refineSRBp(const Options* opt, const bam_hdr_t* hdr, const char* liteChrSeq, const char* largeChrSeq){
    if((int32_t)mConsensus.size() < 2 * opt->filterOpt->mMinFlankSize) return false;
    // Get reference slice
    BreakPoint bp = BreakPoint(this, hdr);
    if(mSVT >= 5) bp = BreakPoint(this, hdr, 10 * opt->libInfo->mReadLen);
    if(liteChrSeq || largeChrSeq) mSVRef = bp.getSVRef(liteChrSeq, largeChrSeq);
    else mergeRef();
    // SR consensus to mSVRef alignment
    Matrix2D<char>* alnResult = new Matrix2D<char>();
    if(!consensusRefAlign(alnResult)){
        delete alnResult;
        return false;
    }
    // Check breakpoint
    AlignDescriptor ad;
    if(!findSplit(alnResult, ad)){
#ifdef DEBUG
        if(opt->debug & DEBUG_FCALL){
            std::cout << "debug_find_split_fail_SVID: " << mID << std::endl;
            std::cout << *alnResult << std::endl;
            std::cout << ad << std::endl;
        }
#endif
        delete alnResult;
        return false;
    }
    delete alnResult;
    // Get the start and end of the SV
    int32_t finalGapStart = 0;
    int32_t finalGapEnd = 0;
    if(!coordTransform(bp, ad, finalGapStart, finalGapEnd)) return false;
    if(mSVT < 4 && mSVStart >= mSVEnd) return false;
    // Precise SV from SR and SVRef split alignment found
    mPrecise = true;
    mSVStart = std::max(0, finalGapStart);
    mSVEnd = std::min(bp.mSVEndEnd, finalGapEnd);
    mSRAlignQuality = ad.mPercID;
    mAlnInsLen = ad.mCSEnd - ad.mCSStart - 1;
    mHomLen = std::max(0, ad.mHomLeft + ad.mHomRight - 2);
    int32_t ciWiggle = std::max(ad.mHomLeft, ad.mHomRight);
    mCiPosLow = -ciWiggle;
    mCiPosHigh = ciWiggle;
    mCiEndLow = -ciWiggle;
    mCiEndHigh = ciWiggle;
    // Store gap coordinates
    mGapCoord[0] = ad.mCSStart;
    mGapCoord[1] = ad.mCSEnd;
    mGapCoord[2] = ad.mRefStart;
    mGapCoord[3] = ad.mRefEnd;
    mProbeBegC = mConsensus.substr(std::max(0, ad.mCSStart - opt->filterOpt->mMinFlankSize), 2 * opt->filterOpt->mMinFlankSize);
    mProbeEndC = mConsensus.substr(std::max(0, ad.mCSEnd - opt->filterOpt->mMinFlankSize), 2 * opt->filterOpt->mMinFlankSize);
    // Get probe sequences used for allele fraction computation
    if(mBpInsSeq.length() > 0){// Put inserted sequence after breakpoint back if possible
        mConsensus = mConsensus.substr(0, ad.mCSStart) + mBpInsSeq + mConsensus.substr(ad.mCSStart);
        ad.mCSEnd += mBpInsSeq.length();
        mProbeBegA = mConsensus.substr(std::max(0, ad.mCSStart - opt->filterOpt->mMinFlankSize), 2 * opt->filterOpt->mMinFlankSize);
        mProbeEndA = mConsensus.substr(std::max(0, ad.mCSEnd - opt->filterOpt->mMinFlankSize), 2 * opt->filterOpt->mMinFlankSize);
    }
    // Get inserted sequence if it's an insertion event
    if(mSVT == 4) mInsSeq = mConsensus.substr(ad.mCSStart, ad.mCSEnd - ad.mCSStart - 1);
    // Consensus of SRs and reference aligned successfully
    return true;
}

void mergeSRSVs(SVSet& sr, SVSet& msr, Options* opt){
    // repeat region filter and bp refining
    util::loginfo("Beg online BWA realignment");
    std::vector<std::future<int>> alnret(sr.size());
    for(uint32_t i = 0; i < sr.size(); ++i){
        if(sr[i]->mSVT == 4) continue;
        alnret[i] = opt->pool->enqueue(&RealnFilter::validCCSeq, opt->realnf, std::ref(sr[i]->mConsensus), std::ref(sr[i]->mNameChr1), std::ref(sr[i]->mSVStart), std::ref(sr[i]->mNameChr2), std::ref(sr[i]->mSVEnd), sr[i]->mGapCoord[0], sr[i]->mBpInsSeq.size());
    }
    for(uint32_t i = 0; i < alnret.size(); ++i){
        if(sr[i]->mSVT != 4){
            sr[i]->mRealnRet = alnret[i].get();
        }
    }
    util::loginfo("End online BWA realignment");
#ifdef DEBUG
    if(opt->debug & DEBUG_FREAN){
        std::cout << "\ndebug_realign_failed_sv_info:" << std::endl;
        for(uint32_t i = 0; i < sr.size(); ++i){
            if(sr[i]->mRealnRet < 0 || sr[i]->mRealnRet > 4){
                std::cout << sr[i] << std::endl;
            }
        }
    }
#endif
    // sort 
    std::sort(sr.begin(), sr.end(), SortSVOne());
    for(uint32_t i = 0; i < sr.size(); ++i) sr[i]->mID = i;
    util::loginfo("Beg merging SR supported SVs, raw " + std::to_string(sr.size()));
    int32_t maxCI = 0;
    for(uint32_t i = 0; i < sr.size(); ++i){
        maxCI = std::max(maxCI, std::max(std::abs(sr[i]->mCiPosLow), std::abs(sr[i]->mCiPosHigh)));
        maxCI = std::max(maxCI, std::max(std::abs(sr[i]->mCiEndLow), std::abs(sr[i]->mCiEndHigh)));
    }
    maxCI = std::max(opt->filterOpt->mMaxReadSep, maxCI);
#ifdef DEBUG
    if(opt->debug & DEBUG_FCALL){
        std::cout << "debug_maxCI_used_in_merge_SRSVS: " << maxCI << std::endl;
    }
#endif
    int32_t totSV = sr.size();
    for(int32_t i = 0; i < totSV; ++i){
        if(sr[i]->mMerged) continue;
        if(sr[i]->mSVT == 4 && sr[i]->mInsSeq.length() > 0){ // Keep all insertions
            msr.push_back(sr[i]);
            sr[i]->mMerged = true;
        }
        if(sr[i]->mSRSupport == 0 || sr[i]->mSRAlignQuality == 0){
            sr[i]->mMerged = true;
            continue; // SR assembly failed
        }
        for(int32_t j = i - 1; j >= 0; --j){
            if(sr[j]->mMerged) continue;
            if(sr[i]->mSVT != sr[j]->mSVT || sr[i]->mChr1 != sr[j]->mChr1 || sr[i]->mChr2 != sr[j]->mChr2) break;
            if(std::abs(sr[j]->mSVStart - sr[i]->mSVStart) > maxCI) break;
            // Test whether breakpoints within SR condidence interval
            if(sr[i]->mSVStart >= sr[j]->mSVStart - maxCI && sr[i]->mSVStart <= sr[j]->mSVStart + maxCI &&
               sr[i]->mSVEnd >= sr[j]->mSVEnd - maxCI && sr[i]->mSVEnd <= sr[j]->mSVEnd + maxCI){
                if(sr[i]->mSRSupport < sr[j]->mSRSupport || (i < j && sr[i]->mSRSupport == sr[j]->mSRSupport)) sr[i]->mMerged = true;
            }
        }
        for(int32_t j = i + 1; j < totSV; ++j){
            if(sr[j]->mMerged) continue;
            if(sr[i]->mSVT != sr[j]->mSVT || sr[i]->mChr1 != sr[j]->mChr1 || sr[i]->mChr2 != sr[j]->mChr2) break;
            if(std::abs(sr[j]->mSVStart - sr[i]->mSVStart) > maxCI) break;
            // Test whether breakpoints within SR condidence interval
            if(sr[i]->mSVStart >= sr[j]->mSVStart - maxCI && sr[i]->mSVStart <= sr[j]->mSVStart + maxCI &&
               sr[i]->mSVEnd >= sr[j]->mSVEnd - maxCI && sr[i]->mSVEnd <= sr[j]->mSVEnd + maxCI){
                if(sr[i]->mSRSupport < sr[j]->mSRSupport || (i < j && sr[i]->mSRSupport == sr[j]->mSRSupport)) sr[i]->mMerged = true;
            }
        }
    }
    std::copy_if(sr.begin(), sr.end(), std::back_inserter(msr), [&](const SVRecord* sv){return !sv->mMerged;});
    util::loginfo("End merging SR supported SVs, got " + std::to_string(msr.size()));
#ifdef DEBUG
    if(opt->debug & DEBUG_FCALL){
        std::cout << "debug_Merged_SR_SV_IDs:";
        for(uint32_t i = 0; i < sr.size(); ++i){
            if(sr[i]->mMerged){
                std::cout << i << "\t";
            }
        }
        std::cout << std::endl;
    }
#endif
    // clean merged ones and set unmerge to NULL in case freed
    for(auto& e: sr){
        if(e->mMerged) delete e;
        e = NULL;
    }
    sr.clear(); sr.shrink_to_fit();
}

void mergeDPSVs(SVSet& dp, SVSet& mdp, Options* opt){
    // sort
    std::sort(dp.begin(), dp.end(), SortSVOne());
    for(uint32_t i = 0; i < dp.size(); ++i) dp[i]->mID = i;
    // filter invalid dp svs firstly
    for(uint32_t dpi = 0; dpi < dp.size(); ++dpi){
        if(dp[dpi]->mSVStart < 0 || dp[dpi]->mSVEnd < 0 ||
           dp[dpi]->mSVEnd >= (int32_t)opt->bamheader->target_len[dp[dpi]->mChr2] ||
           dp[dpi]->mSVStart >= (int32_t)opt->bamheader->target_len[dp[dpi]->mChr1]){
            dp[dpi]->mMerged = true;
        }
    }
    // then do merge
    util::loginfo("Beg merging DP supported SVs, raw " + std::to_string(dp.size()));
    // index dpsvs
    int32_t totSV = dp.size();
    for(int32_t i = 0; i < totSV; ++i){
        if(dp[i]->mMerged) continue;
        for(int32_t j = i - 1; j >= 0; --j){
            if(dp[j]->mMerged) continue;
            if(dp[i]->mSVT != dp[j]->mSVT || dp[i]->mChr1 != dp[j]->mChr1 || dp[i]->mChr2 != dp[j]->mChr2) break;
            if(std::abs(dp[j]->mSVStart - dp[i]->mSVStart) > opt->libInfo->mMaxNormalISize) break;
            // Test whether breakpoints within SR condidence interval
            if((std::abs(dp[i]->mSVStart - dp[j]->mSVStart) < opt->libInfo->mMaxNormalISize) &&
               (std::abs(dp[i]->mSVEnd - dp[j]->mSVEnd) < opt->libInfo->mMaxNormalISize)){
                if(dp[i]->mPESupport < dp[j]->mPESupport || (i < j && dp[i]->mPESupport == dp[j]->mPESupport)) dp[i]->mMerged = true;
            }
        }
        for(int32_t j = i + 1; j < totSV; ++j){
            if(dp[j]->mMerged) continue;
            if(dp[i]->mSVT != dp[j]->mSVT || dp[i]->mChr1 != dp[j]->mChr1 || dp[i]->mChr2 != dp[j]->mChr2) break;
            if(std::abs(dp[j]->mSVStart - dp[i]->mSVStart) > opt->libInfo->mMaxNormalISize) break;
            // Test whether breakpoints within SR condidence interval
            if((std::abs(dp[i]->mSVStart - dp[j]->mSVStart) < opt->libInfo->mMaxNormalISize) &&
               (std::abs(dp[i]->mSVEnd - dp[j]->mSVEnd) < opt->libInfo->mMaxNormalISize)){
                if(dp[i]->mPESupport < dp[j]->mPESupport || (i < j && dp[i]->mPESupport == dp[j]->mPESupport)) dp[i]->mMerged = true;
            }
        }
    }
    std::copy_if(dp.begin(), dp.end(), std::back_inserter(mdp), [&](const SVRecord* sv){return !sv->mMerged;});
    util::loginfo("End merging DP supported SVs, got " + std::to_string(mdp.size()));
#ifdef DEBUG
    if(opt->debug & DEBUG_FCALL){
        std::cout << "debug_Merged_DP_SV_IDs:";
        for(uint32_t i = 0; i < dp.size(); ++i){
            if(dp[i]->mMerged){
                std::cout << "\t" << i;
            }
        }
        std::cout << std::endl;
    }
#endif
    // clean merged ones and set unmerge to NULL in case freed
    for(auto& e: dp){
        if(e->mMerged) delete e;
        e = NULL;
    }
    dp.clear(); dp.shrink_to_fit();
}

void mergeAndSortSVSet(SVSet& sr, SVSet& dp, SVSet& svs, Options* opt){
    // Merge SR SVSet
    mergeSRSVs(sr, svs, opt);
    // Merge DP SVset
    SVSet pe;
    mergeDPSVs(dp, pe, opt);
    // index dpsvs
    std::vector<std::vector<std::pair<int32_t, int32_t>>> chridx(9);
    for(uint32_t i = 0; i < 9; ++i) chridx[i].resize(opt->contigNum, {-1, -1});
    for(int32_t i = 0; i < (int32_t)pe.size(); ++i){
        if(chridx[pe[i]->mSVT][pe[i]->mChr1].first < 0) chridx[pe[i]->mSVT][pe[i]->mChr1].first = i;
        chridx[pe[i]->mSVT][pe[i]->mChr1].second = i;
    }
    // Augment SR SVs with PE records
    util::loginfo("Beg augmenting SR supported SV candidates with DPs");
    std::map<int32_t, std::vector<int32_t>> pemset;
    int32_t ttsrs = svs.size();
    for(int32_t i = 0; i < ttsrs; ++i){
        if(chridx[svs[i]->mSVT][svs[i]->mChr1].first < 0) continue;
        for(int32_t j = chridx[svs[i]->mSVT][svs[i]->mChr1].first; j <= chridx[svs[i]->mSVT][svs[i]->mChr1].second; ++j){
            if(svs[i]->mChr2 > pe[j]->mChr2) continue;
            if(svs[i]->mChr2 < pe[j]->mChr2) break;
            // Test whether breakpoint is within PE confidence interval
            if(svs[i]->mSVStart >= pe[j]->mSVStart - std::max(opt->libInfo->mMaxNormalISize, pe[j]->mCiPosHigh) && 
               svs[i]->mSVStart <= pe[j]->mSVStart + std::max(opt->libInfo->mMaxNormalISize, pe[j]->mCiPosHigh) &&
               svs[i]->mSVEnd >= pe[j]->mSVEnd - std::max(opt->libInfo->mMaxNormalISize, pe[j]->mCiEndHigh) && 
               svs[i]->mSVEnd <= pe[j]->mSVEnd + std::max(opt->libInfo->mMaxNormalISize, pe[j]->mCiEndHigh)){
                auto iter = pemset.find(j);
                if(iter == pemset.end()) pemset[j] = std::vector<int32_t>(1, i);
                else iter->second.push_back(i);
            }
        }
    }
    for(auto iter = pemset.begin(); iter != pemset.end(); ++iter){
        int32_t pmid = iter->second[0];
        if(iter->second.size() > 1){
            std::vector<int32_t> ndps;
            for(auto i: iter->second){
                if(svs[i]->mRealnRet >= 0 && svs[i]->mRealnRet <= opt->fuseOpt->mWhiteFilter.mMaxRepHit){
                    ndps.push_back(i);
                }
            }
            if(ndps.size() == 1){
                pmid = ndps[0];
            }else if(ndps.size() > 2){
                pmid = ndps[0];
                for(uint32_t j = 1; j < ndps.size(); ++j){
                    if(svs[ndps[j]]->mSRSupport > svs[pmid]->mSRSupport){
                        pmid = ndps[j];
                    }
                }
            }else if(ndps.empty()){
                for(uint32_t j = 1; j < iter->second.size(); ++j){
                    if(svs[iter->second[j]]->mSRSupport > svs[pmid]->mSRSupport){
                        pmid = iter->second[j];
                    }
                }
            }
        }
        svs[pmid]->mPESupport = pe[iter->first]->mPESupport;
        svs[pmid]->mPEMapQuality = pe[iter->first]->mPEMapQuality;
        pe[iter->first]->mMerged = true;
    }
    util::loginfo("End augmenting SR supported SV candidates with DPs");
    // Append DP supported only SV 
    std::copy_if(pe.begin(), pe.end(), std::back_inserter(svs), [&](const SVRecord* sv){return !sv->mMerged;});
    // clean merged ones and set unmerge to NULL in case freed
    for(auto& e: pe){
        if(e->mMerged) delete e;
        e = NULL;
    }
    pe.clear(); pe.shrink_to_fit();
}

void getDPSVRef(SVSet& pe, Options* opt){
    // Open file handler
    faidx_t* fai = fai_load(opt->alnref.c_str());
    // get SVRef on same chr
    for(auto& e: pe){
        if(e->mPrecise) continue;
        e->mNameChr1 = opt->bamheader->target_name[e->mChr1];
        e->mNameChr2 = opt->bamheader->target_name[e->mChr2];
        int32_t chrLen = -1;
        char* chrSeq = faidx_fetch_seq(fai, opt->bamheader->target_name[e->mChr1],  e->mSVStart - 1, e->mSVStart - 1, &chrLen);
        e->mSVRef = std::string(1, std::toupper(chrSeq[0]));
        free(chrSeq);
    }
    fai_destroy(fai);
}
