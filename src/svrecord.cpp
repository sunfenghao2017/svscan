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

void mergeSRSVs(SVSet* sr, SVSet* msr, Options* opt){
    // repeat region filter and bp refining
    util::loginfo("Beg online BWA realignment");
    std::vector<std::future<int>> alnret(sr->size());
    for(uint32_t i = 0; i < sr->size(); ++i){
        if(sr->at(i)->mSVT == 4) continue;
        alnret[i] = opt->pool->enqueue(&RealnFilter::validCCSeq, opt->realnf, std::ref(sr->at(i)->mConsensus), std::ref(sr->at(i)->mNameChr1), std::ref(sr->at(i)->mSVStart), std::ref(sr->at(i)->mNameChr2), std::ref(sr->at(i)->mSVEnd), sr->at(i)->mGapCoord[0], sr->at(i)->mBpInsSeq.size());
    }
    for(uint32_t i = 0; i < alnret.size(); ++i){
        if(sr->at(i)->mSVT != 4){
            sr->at(i)->mRealnRet = alnret[i].get();
        }
    }
    util::loginfo("End online BWA realignment");
#ifdef DEBUG
    if(opt->debug & DEBUG_FREAN){
        std::cout << "\ndebug_realign_failed_sv_info:" << std::endl;
        for(uint32_t i = 0; i < sr->size(); ++i){
            if(sr->at(i)->mRealnRet < 0 || sr->at(i)->mRealnRet > 4){
                std::cout << sr[i] << std::endl;
            }
        }
    }
#endif
    // sort 
    std::sort(sr->begin(), sr->end(), SortSVOne());
    for(uint32_t i = 0; i < sr->size(); ++i) sr->at(i)->mID = i;
    util::loginfo("Beg merging SR supported SVs, raw " + std::to_string(sr->size()));
    int32_t maxCI = 0;
    for(uint32_t i = 0; i < sr->size(); ++i){
        maxCI = std::max(maxCI, std::max(std::abs(sr->at(i)->mCiPosLow), std::abs(sr->at(i)->mCiPosHigh)));
        maxCI = std::max(maxCI, std::max(std::abs(sr->at(i)->mCiEndLow), std::abs(sr->at(i)->mCiEndHigh)));
    }
    maxCI = std::max(opt->filterOpt->mMaxReadSep, maxCI);
#ifdef DEBUG
    if(opt->debug & DEBUG_FCALL){
        std::cout << "debug_maxCI_used_in_merge_SRSVS: " << maxCI << std::endl;
    }
#endif
    int32_t totSV = sr->size();
    for(int32_t i = 0; i < totSV; ++i){
        if(sr->at(i)->mMerged) continue;
        if(sr->at(i)->mSVT == 4 && sr->at(i)->mInsSeq.length() > 0){ // Keep all insertions
            msr->push_back(sr->at(i));
            sr->at(i)->mMerged = true;
        }
        if(sr->at(i)->mSRSupport == 0 || sr->at(i)->mSRAlignQuality == 0){
            sr->at(i)->mMerged = true;
            continue; // SR assembly failed
        }
        for(int32_t j = i - 1; j >= 0; --j){
            if(sr->at(j)->mMerged) continue;
            if(sr->at(i)->mSVT != sr->at(j)->mSVT || sr->at(i)->mChr1 != sr->at(j)->mChr1 || sr->at(i)->mChr2 != sr->at(j)->mChr2) break;
            if(std::abs(sr->at(j)->mSVStart - sr->at(i)->mSVStart) > maxCI) break;
            // Test whether breakpoints within SR condidence interval
            if(sr->at(i)->mSVStart >= sr->at(j)->mSVStart - maxCI && sr->at(i)->mSVStart <= sr->at(j)->mSVStart + maxCI &&
               sr->at(i)->mSVEnd >= sr->at(j)->mSVEnd - maxCI && sr->at(i)->mSVEnd <= sr->at(j)->mSVEnd + maxCI){
                if(sr->at(i)->mSRSupport < sr->at(j)->mSRSupport || (i < j && sr->at(i)->mSRSupport == sr->at(j)->mSRSupport)) sr->at(i)->mMerged = true;
            }
        }
        for(int32_t j = i + 1; j < totSV; ++j){
            if(sr->at(j)->mMerged) continue;
            if(sr->at(i)->mSVT != sr->at(j)->mSVT || sr->at(i)->mChr1 != sr->at(j)->mChr1 || sr->at(i)->mChr2 != sr->at(j)->mChr2) break;
            if(std::abs(sr->at(j)->mSVStart - sr->at(i)->mSVStart) > maxCI) break;
            // Test whether breakpoints within SR condidence interval
            if(sr->at(i)->mSVStart >= sr->at(j)->mSVStart - maxCI && sr->at(i)->mSVStart <= sr->at(j)->mSVStart + maxCI &&
               sr->at(i)->mSVEnd >= sr->at(j)->mSVEnd - maxCI && sr->at(i)->mSVEnd <= sr->at(j)->mSVEnd + maxCI){
                if(sr->at(i)->mSRSupport < sr->at(j)->mSRSupport || (i < j && sr->at(i)->mSRSupport == sr->at(j)->mSRSupport)) sr->at(i)->mMerged = true;
            }
        }
    }
    for(uint32_t i = 0; i < sr->size(); ++i){
        if(sr->at(i)->mMerged){
            delete sr->at(i);
            sr->at(i) = NULL;
        }else{
            msr->push_back(sr->at(i));
            sr->at(i) = NULL;
        }
    }
    util::loginfo("End merging SR supported SVs, got " + std::to_string(msr->size()));
#ifdef DEBUG
    if(opt->debug & DEBUG_FCALL){
        std::cout << "debug_Merged_SR_SV_IDs:";
        for(uint32_t i = 0; i < sr.size(); ++i){
            if(sr->at(i)->mMerged){
                std::cout << i << "\t";
            }
        }
        std::cout << std::endl;
    }
#endif
    sr->clear(); sr->shrink_to_fit();
}

void mergeDPSVs(SVSet* dp, SVSet* mdp, Options* opt){
    // sort
    std::sort(dp->begin(), dp->end(), SortSVOne());
    for(uint32_t i = 0; i < dp->size(); ++i) dp->at(i)->mID = i;
    // filter invalid dp svs firstly
    for(uint32_t dpi = 0; dpi < dp->size(); ++dpi){
        if(dp->at(dpi)->mSVStart < 0 || dp->at(dpi)->mSVEnd < 0 ||
           dp->at(dpi)->mSVEnd >= (int32_t)opt->bamheader->target_len[dp->at(dpi)->mChr2] ||
           dp->at(dpi)->mSVStart >= (int32_t)opt->bamheader->target_len[dp->at(dpi)->mChr1]){
            dp->at(dpi)->mMerged = true;
        }
    }
    // then do merge
    util::loginfo("Beg merging DP supported SVs, raw " + std::to_string(dp->size()));
    int32_t totSV = dp->size();
    for(int32_t i = 0; i < totSV; ++i){
        if(dp->at(i)->mMerged) continue;
        for(int32_t j = i - 1; j >= 0; --j){
            if(dp->at(j)->mMerged) continue;
            if(dp->at(i)->mSVT != dp->at(j)->mSVT || dp->at(i)->mChr1 != dp->at(j)->mChr1 || dp->at(i)->mChr2 != dp->at(j)->mChr2) break;
            if(std::abs(dp->at(j)->mSVStart - dp->at(i)->mSVStart) > opt->libInfo->mMaxNormalISize) break;
            // Test whether breakpoints within SR condidence interval
            if((std::abs(dp->at(i)->mSVStart - dp->at(j)->mSVStart) < opt->libInfo->mMaxNormalISize) &&
               (std::abs(dp->at(i)->mSVEnd - dp->at(j)->mSVEnd) < opt->libInfo->mMaxNormalISize)){
                if(dp->at(i)->mPESupport < dp->at(j)->mPESupport || (i < j && dp->at(i)->mPESupport == dp->at(j)->mPESupport)) dp->at(i)->mMerged = true;
            }
        }
        for(int32_t j = i + 1; j <= totSV; ++j){
            if(dp->at(j)->mMerged) continue;
            if(dp->at(i)->mSVT != dp->at(j)->mSVT || dp->at(i)->mChr1 != dp->at(j)->mChr1 || dp->at(i)->mChr2 != dp->at(j)->mChr2) break;
            if(std::abs(dp->at(j)->mSVStart - dp->at(i)->mSVStart) > opt->libInfo->mMaxNormalISize) break;
            // Test whether breakpoints within SR condidence interval
            if((std::abs(dp->at(i)->mSVStart - dp->at(j)->mSVStart) < opt->libInfo->mMaxNormalISize) &&
               (std::abs(dp->at(i)->mSVEnd - dp->at(j)->mSVEnd) < opt->libInfo->mMaxNormalISize)){
                if(dp->at(i)->mPESupport < dp->at(j)->mPESupport || (i < j && dp->at(i)->mPESupport == dp->at(j)->mPESupport)) dp->at(i)->mMerged = true;
            }
        }
    }
    for(uint32_t i = 0; i < dp->size(); ++i){
        if(dp->at(i)->mMerged){
            delete dp->at(i);
            dp->at(i) = NULL;
        }else{
            mdp->push_back(dp->at(i));
            dp->at(i) = NULL;
        }
    }
    util::loginfo("End merging DP supported SVs, got " + std::to_string(mdp->size()));
#ifdef DEBUG
    if(opt->debug & DEBUG_FCALL){
        std::cout << "debug_Merged_DP_SV_IDs:";
        for(uint32_t i = 0; i < dp.size(); ++i){
            if(dp->at(i)->mMerged){
                std::cout << "\t" << i;
            }
        }
        std::cout << std::endl;
    }
#endif
    dp->clear(); dp->shrink_to_fit();
}

void mergeAndSortSVSet(SVSet* sr, SVSet* dp, SVSet* svs, Options* opt){
    // Merge SR SVSet
    mergeSRSVs(sr, svs, opt);
    // Merge DP SVset
    SVSet pe;
    mergeDPSVs(dp, &pe, opt);
    // index dpsvs
    std::vector<std::vector<std::pair<int32_t, int32_t>>> chridx(9);
    for(uint32_t i = 0; i < 9; ++i){
        chridx[i].resize(opt->contigNum, {-1, -1});
    }
    for(int32_t i = 0; i < (int32_t)dp->size(); ++i){
        if(chridx[dp->at(i)->mSVT][dp->at(i)->mChr1].first < 0) chridx[dp->at(i)->mSVT][dp->at(i)->mChr1].first = i;
        chridx[dp->at(i)->mSVT][dp->at(i)->mChr1].second = i;
    }
    // Augment SR SVs with PE records
    util::loginfo("Beg augmenting SR supported SV candidates with DPs");
    std::map<int32_t, std::vector<int32_t>> pemset;
    int32_t ttsrs = svs->size();
    for(int32_t i = 0; i < ttsrs; ++i){
        if(chridx[svs->at(i)->mSVT][svs->at(i)->mChr1].first < 0) continue;
        for(int32_t j = chridx[svs->at(i)->mSVT][svs->at(i)->mChr1].first; j <= chridx[svs->at(i)->mSVT][svs->at(i)->mChr1].second; ++j){
            if(pe[j]->mChr2 > svs->at(i)->mChr2) break;
            if(pe[j]->mChr2 < svs->at(i)->mChr2) continue;
            // Test whether breakpoint is within PE confidence interval
            if(svs->at(i)->mSVStart >= pe[j]->mSVStart - std::max(opt->libInfo->mMaxNormalISize, pe[j]->mCiPosHigh) && 
               svs->at(i)->mSVStart <= pe[j]->mSVStart + std::max(opt->libInfo->mMaxNormalISize, pe[j]->mCiPosHigh) &&
               svs->at(i)->mSVEnd >= pe[j]->mSVEnd - std::max(opt->libInfo->mMaxNormalISize, pe[j]->mCiEndHigh) && 
               svs->at(i)->mSVEnd <= pe[j]->mSVEnd + std::max(opt->libInfo->mMaxNormalISize, pe[j]->mCiEndHigh)){
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
                if(svs->at(i)->mRealnRet >= 0 && svs->at(i)->mRealnRet <= opt->fuseOpt->mWhiteFilter.mMaxRepHit){
                    ndps.push_back(i);
                }
            }
            if(ndps.size() == 1){
                pmid = ndps[0];
            }else if(ndps.size() > 2){
                pmid = ndps[0];
                for(uint32_t j = 1; j < ndps.size(); ++j){
                    if(svs->at(ndps[j])->mSRSupport > svs->at(pmid)->mSRSupport){
                        pmid = ndps[j];
                    }
                }
            }else if(ndps.empty()){
                for(uint32_t j = 1; j < iter->second.size(); ++j){
                    if(svs->at(iter->second[j])->mSRSupport > svs->at(pmid)->mSRSupport){
                        pmid = iter->second[j];
                    }
                }
            }
        }
        svs->at(pmid)->mPESupport = pe[iter->first]->mPESupport;
        svs->at(pmid)->mPEMapQuality = pe[iter->first]->mPEMapQuality;
        pe[iter->first]->mMerged = true;
    }
    util::loginfo("End augmenting SR supported SV candidates with DPs");
    // Append DP supported only SV 
    for(uint32_t i = 0; i < pe.size(); ++i){
        if(pe[i]->mMerged){
            delete pe[i];
            pe[i] = NULL;
        }else{
            svs->push_back(pe[i]);
            pe[i] = NULL;
        }
    }
    pe.resize(0); pe.shrink_to_fit();
}

void getDPSVRef(SVSet* pe, Options* opt){
    // Open file handler
    faidx_t* fai = fai_load(opt->alnref.c_str());
    // get SVRef on same chr
    for(uint32_t i = 0; i < pe->size(); ++i){
        if(pe->at(i)->mPrecise) continue;
        pe->at(i)->mNameChr1 = opt->bamheader->target_name[pe->at(i)->mChr1];
        pe->at(i)->mNameChr2 = opt->bamheader->target_name[pe->at(i)->mChr2];
        int32_t chrLen = -1;
        char* chrSeq = faidx_fetch_seq(fai, opt->bamheader->target_name[pe->at(i)->mChr1],  pe->at(i)->mSVStart - 1, pe->at(i)->mSVStart - 1, &chrLen);
        pe->at(i)->mSVRef = std::string(1, std::toupper(chrSeq[0]));
        free(chrSeq);
    }
    fai_destroy(fai);
}
