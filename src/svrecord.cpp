#include "svrecord.h"
#include "breakpoint.h"

template<typename TBreakPoint>
bool SVRecord::coordTransform(TBreakPoint& bp, AlignDescriptor& ad, int32_t& finalGapStart, int32_t& finalGapEnd){
    if(mSVT == 0 || mSVT == 5){// 5to5
        int32_t annealed = bp.mSVStartEnd - bp.mSVStartBeg;
        if(ad.mRefStart >= annealed || (ad.mRefEnd < annealed)) return false;
        finalGapStart = bp.mSVStartBeg + ad.mRefStart;
        finalGapEnd = bp.mSVEndBeg + mSVRef.size() - ad.mRefEnd + 1;
        return true;
    }
    if(mSVT == 1 || mSVT == 6){// 3to3
        int32_t annealed = bp.mSVStartEnd - bp.mSVStartBeg;
        if(ad.mRefStart >= annealed || (ad.mRefEnd < annealed)) return false;
        finalGapStart = bp.mSVStartBeg + (annealed - ad.mRefStart) - 1;
        finalGapEnd = bp.mSVEndBeg + (ad.mRefEnd - annealed) - 2;
        return true;
    }
    if(mSVT == 2 || mSVT == 7){// 5to3
        int32_t annealed = bp.mSVStartEnd - bp.mSVStartBeg;
        if(ad.mRefStart >= annealed || (ad.mRefEnd < annealed)) return false;
        finalGapStart = bp.mSVStartBeg + ad.mRefStart;
        finalGapEnd = bp.mSVEndBeg + (ad.mRefEnd - annealed) - 2;
        return true;
    }
    if(mSVT == 3 || mSVT == 8){// 3to5
        int32_t annealed = bp.mSVEndEnd - bp.mSVEndBeg;
        if(ad.mRefStart >= annealed || (ad.mRefEnd < annealed)) return false;
        finalGapStart = bp.mSVStartBeg + (ad.mRefEnd - annealed) - 2;
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

bool SVRecord::refineSRBp(const Options* opt, const bam_hdr_t* hdr, const char* chr1Seq, const char* chr2Seq){
    if((int32_t)mConsensus.size() < 2 * opt->filterOpt->mMinFlankSize) return false;
    // Get reference slice
    BreakPoint bp = BreakPoint(*this, hdr);
    mSVRef = bp.getSVRef(chr1Seq, chr2Seq);
    // SR consensus to mSVRef alignment
    Matrix2D<char>* alnResult = new Matrix2D<char>();
    if(!consensusRefAlign(alnResult)){
        delete alnResult;
        return false;
    }
    // Check breakpoint
    AlignDescriptor ad;
    if(!findSplit(alnResult, ad)){
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
    mSVStart = finalGapStart;
    mSVEnd = finalGapEnd;
    mSRAlignQuality = ad.mPercID;
    mInsLen = ad.mCSEnd - ad.mCSStart - 1;
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
    // Get probe sequences used for allele fraction computation
    mProbeBegC = mConsensus.substr(std::max(0, ad.mCSStart - opt->filterOpt->mMinFlankSize), 2 * opt->filterOpt->mMinFlankSize + 1);
    mProbeBegR = mSVRef.substr(std::max(0, ad.mRefStart -opt->filterOpt->mMinFlankSize), std::min(2 * opt->filterOpt->mMinFlankSize + 1, bp.mSVStartEnd - bp.mSVStartBeg));
    mProbeEndC = mConsensus.substr(std::max(0, ad.mCSEnd - opt->filterOpt->mMinFlankSize), 2 * opt->filterOpt->mMinFlankSize + 1);
    mProbeEndR = mSVRef.substr(std::max(0, ad.mRefEnd - std::min(bp.mSVStartEnd - bp.mSVStartBeg, opt->filterOpt->mMinFlankSize)), 2 * opt->filterOpt->mMinFlankSize + 1);
    // Get inserted sequence if it's an insertion event
    if(mSVT == 4) mInsSeq = mConsensus.substr(ad.mCSStart, ad.mCSEnd - ad.mCSStart - 1);
    // Consensus of SRs and reference aligned successfully
    return true;
}

void mergeAndSortSVSet(SVSet& sr, SVSet& pe, int32_t peMergeSearchWindow, int32_t srDupSearchWindow){
    // Sort PE records for look-up
    util::loginfo("Start sorting DP supported SV candidates");
    std::sort(pe.begin(), pe.end());
    util::loginfo("Finish sorting DP supported SV candidates");
    // Sort SR records for look-up
    util::loginfo("Start sorting SR supported SV candidates");
    std::sort(sr.begin(), sr.end());
    util::loginfo("Finish sorting SR supported SV candidates");
    // Augment PE SVs and append missing SR SVs
    util::loginfo("Start search DP supports for each SR supported SV candidates");
    for(int32_t i = 0; i < (int32_t)sr.size(); ++i){
        util::loginfo("Start processing SR supported candidate: " + std::to_string(i));
        if(sr[i].mSRSupport == 0 || sr[i].mSRAlignQuality == 0) continue; // SR assembly failed
        if(sr[i].mSVT == 4){// append insertion anyway
            pe.push_back(sr[i]);
            continue;
        }
        // Precise duplicates
        bool svExists = false;
        SVRecord searchSV;
        searchSV.mChr1 = sr[i].mChr1;
        searchSV.mSVStart = sr[i].mSVStart - peMergeSearchWindow;
        searchSV.mSVEnd = sr[i].mSVEnd;
        auto itOther = std::lower_bound(pe.begin(), pe.end(), searchSV);
        for(;itOther != pe.end() && (std::abs(itOther->mSVStart - sr[i].mSVStart) < peMergeSearchWindow); ++itOther){
            if(itOther->mSVT != sr[i].mSVT || itOther->mPrecise) continue;
            if(sr[i].mChr1 != itOther->mChr1 || sr[i].mChr2 != itOther->mChr2) continue;
            // Test whether breakpoint is within PE confidence interval
            if((itOther->mSVStart + itOther->mCiPosLow) < sr[i].mSVStart && (itOther->mSVStart + itOther->mCiPosHigh) > sr[i].mSVStart &&
               (itOther->mSVEnd + itOther->mCiEndLow) < sr[i].mSVEnd && (itOther->mSVEnd + itOther->mCiEndHigh) > sr[i].mSVEnd){
                svExists = true;
                // Augment PE record
                itOther->mSVStart = sr[i].mSVStart;
                itOther->mSVEnd = sr[i].mSVEnd;
                itOther->mCiPosLow = sr[i].mCiPosLow;
                itOther->mCiPosHigh = sr[i].mCiPosHigh;
                itOther->mCiEndLow = sr[i].mCiEndLow;
                itOther->mCiEndHigh = sr[i].mCiEndHigh;
                itOther->mSRSupport = sr[i].mSRSupport;
                itOther->mSRMapQuality = sr[i].mSRMapQuality;
                itOther->mInsLen = sr[i].mInsLen;
                itOther->mHomLen = sr[i].mHomLen;
                itOther->mSRAlignQuality = sr[i].mSRAlignQuality;
                itOther->mPrecise = true;
                itOther->mConsensus = sr[i].mConsensus;
                itOther->mSVRef = sr[i].mSVRef;
                itOther->mProbeBegC = sr[i].mProbeBegC;
                itOther->mProbeBegR = sr[i].mProbeBegR;
                itOther->mProbeEndC = sr[i].mProbeEndC;
                itOther->mProbeEndR = sr[i].mProbeEndR;
                itOther->mNameChr1 = sr[i].mNameChr1;
                itOther->mNameChr2 = sr[i].mNameChr2;
                for(int c = 0; c < 4; ++c) itOther->mGapCoord[c] = sr[i].mGapCoord[c];
            }
        }
        // SR only SV
        if(!svExists){
            // Make sure there is no precise duplicate
            bool preciseDup = false;
            for(int32_t j = i + 1; j < (int32_t)sr.size(); ++j){
                if(std::abs(sr[i].mSVStart - sr[j].mSVStart) > srDupSearchWindow) break;
                if(sr[i].mSVT != sr[j].mSVT) continue;
                if(sr[i].mChr1 != sr[j].mChr1 || sr[i].mChr2 != sr[j].mChr2) continue;
                // Test whether breakpoints within SR condidence interval
                if((sr[j].mSVStart + sr[j].mCiPosLow) <= sr[i].mSVStart && (sr[j].mSVStart + sr[j].mCiPosHigh) >= sr[i].mSVStart &&
                   (sr[j].mSVEnd + sr[j].mCiEndLow) <= sr[i].mSVEnd && (sr[j].mSVEnd + sr[j].mCiPosHigh) >= sr[i].mSVEnd){
                    // Duplicate, keep better call
                    if(sr[i].mSRSupport < sr[j].mSRSupport || (i < j && sr[i].mSRSupport == sr[j].mSRSupport)) preciseDup = true;
                }
            }
            for(int32_t j = i - 1;  j >= 0; --j){
                if(std::abs(sr[i].mSVStart - sr[j].mSVStart) > srDupSearchWindow) break;
                if(sr[i].mSVT != sr[j].mSVT) continue;
                if(sr[i].mChr1 != sr[j].mChr1 || sr[i].mChr2 != sr[j].mChr2) continue;
                // Test whether breakpoints within SR condidence interval
                if((sr[j].mSVStart + sr[j].mCiPosLow) <= sr[i].mSVStart && (sr[j].mSVStart + sr[j].mCiPosHigh) >= sr[i].mSVStart &&
                   (sr[j].mSVEnd + sr[j].mCiEndLow) <= sr[i].mSVEnd && (sr[j].mSVEnd + sr[j].mCiPosHigh) >= sr[i].mSVEnd){
                    // Duplicate, keep better call
                    if(sr[i].mSRSupport < sr[j].mSRSupport || (i < j && sr[i].mSRSupport == sr[j].mSRSupport)) preciseDup = true;
                }

            }
            if(!preciseDup) pe.push_back(sr[i]);
        }
        util::loginfo("Finish processing SR supported candidate: " + std::to_string(i));
    }
    // Re-number SVs
    util::loginfo("Start sorting merged SVs");
    std::sort(pe.begin(), pe.end());
    util::loginfo("Finish sorting merged SVs");
    util::loginfo("Starting reallocating ID of SVs");
    for(uint32_t i = 0; i < pe.size(); ++i) pe[i].mID = i;
    util::loginfo("Finish reallocating ID of SVs");
}

void getDPSVRef(SVSet& pe, Options* opt){
    // Open file handler
    samFile* fp = sam_open(opt->bamfile.c_str(), "r");
    hts_set_fai_filename(fp, opt->genome.c_str());
    bam_hdr_t* h = sam_hdr_read(fp);
    faidx_t* fai = fai_load(opt->genome.c_str());
    // get SVRef on same chr
    for(auto sviter = pe.begin(); sviter != pe.end(); ++sviter){
        if(sviter->mPrecise) continue;
        sviter->mNameChr1 = h->target_name[sviter->mChr1];
        sviter->mNameChr2 = h->target_name[sviter->mChr2];
        int32_t chrLen = -1;
        char* chrSeq = faidx_fetch_seq(fai, h->target_name[sviter->mChr1],  sviter->mSVStart - 1, sviter->mSVStart - 1, &chrLen);
        sviter->mSVRef = std::string(1, std::toupper(chrSeq[0]));
        free(chrSeq);
    }
    sam_close(fp);
    bam_hdr_destroy(h);
    fai_destroy(fai);
}
