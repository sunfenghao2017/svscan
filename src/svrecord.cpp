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
    // Get probe sequences used for allele fraction computation
    mProbeBegR = std::string(chr1Seq + std::max(0, mSVStart - opt->filterOpt->mMinFlankSize), chr1Seq + std::min(bp.mChr1Len, mSVStart + opt->filterOpt->mMinFlankSize));
    mProbeEndR = std::string(chr2Seq + std::max(0, mSVEnd - opt->filterOpt->mMinFlankSize), chr2Seq + std::min(bp.mChr2Len, mSVEnd + opt->filterOpt->mMinFlankSize));
    if(mBpInsSeq.length() > 0){// Put inserted sequence after breakpoint back if possible
        mConsensus = mConsensus.substr(0, ad.mCSStart) + mBpInsSeq + mConsensus.substr(ad.mCSStart);
        ad.mCSEnd += mBpInsSeq.length();
    }
    mProbeBegC = mConsensus.substr(std::max(0, ad.mCSStart - opt->filterOpt->mMinFlankSize), 2 * opt->filterOpt->mMinFlankSize + 1);
    mProbeEndC = mConsensus.substr(std::max(0, ad.mCSEnd - opt->filterOpt->mMinFlankSize), 2 * opt->filterOpt->mMinFlankSize + 1);
    // Get inserted sequence if it's an insertion event
    if(mSVT == 4) mInsSeq = mConsensus.substr(ad.mCSStart, ad.mCSEnd - ad.mCSStart - 1);
    // Consensus of SRs and reference aligned successfully
    return true;
}

void mergeSRSVs(SVSet& sr, SVSet& msr){
    util::loginfo("Starting merge SR supported SVs");
    std::sort(sr.begin(), sr.end());
    for(uint32_t i = 0; i < sr.size(); ++i){
        if(sr[i].mMerged) continue;
        if(sr[i].mSVT == 4 && sr[i].mInsSeq.length() > 0){ // Keep all insertions
            msr.push_back(sr[i]);
            sr[i].mMerged = true;
        }
        if(sr[i].mSRSupport == 0 || sr[i].mSRAlignQuality == 0){
            sr[i].mMerged = true;
            continue; // SR assembly failed
        }
        for(uint32_t j = 0; j < sr.size(); ++j){
            if(i == j || sr[j].mMerged) continue;
            if(sr[i].mSVT != sr[j].mSVT || sr[i].mChr1 != sr[j].mChr1 || sr[i].mChr2 != sr[j].mChr2) continue;
            // Test whether breakpoints within SR condidence interval
            if(sr[i].mSVStart >= sr[j].mSVStart + sr[j].mCiPosLow && sr[i].mSVStart <= sr[j].mSVStart + sr[j].mCiPosHigh &&
               sr[i].mSVEnd >= sr[j].mSVEnd + sr[j].mCiEndLow && sr[i].mSVEnd <= sr[j].mSVEnd + sr[j].mCiEndHigh){
                if(sr[i].mSRSupport < sr[j].mSRSupport || (i < j && sr[i].mSRSupport == sr[j].mSRSupport)) sr[i].mMerged = true;
            }
        }
    }
    std::copy_if(sr.begin(), sr.end(), std::back_inserter(msr), [&](const SVRecord& sv){return !sv.mMerged;});
    util::loginfo(std::to_string(sr.size()) + " SVs merged into " + std::to_string(msr.size()) + " ones");
    sr.clear();
}

void mergeAndSortSVSet(SVSet& sr, SVSet& pe, SVSet& svs){
    // Merge SR SVSet
    mergeSRSVs(sr, svs);
    // Augment SR SVs with PE records
    util::loginfo("Start augmenting SR supported SV candidates with DPs");
    for(uint32_t i = 0; i < svs.size(); ++i){
        for(auto itpe = pe.begin(); itpe != pe.end(); ++itpe){
            if(itpe->mSVT != svs[i].mSVT || svs[i].mChr1 != itpe->mChr1 || svs[i].mChr2 != itpe->mChr2 || itpe->mMerged) continue;
            // Test whether breakpoint is within PE confidence interval
            if(svs[i].mSVStart >= itpe->mSVStart + itpe->mCiPosLow && svs[i].mSVStart <= itpe->mSVStart + itpe->mCiPosHigh &&
               svs[i].mSVEnd >= itpe->mSVEnd + itpe->mCiEndLow && svs[i].mSVEnd <= itpe->mSVEnd + itpe->mCiEndHigh){
                svs[i].mPESupport = itpe->mPESupport;
                svs[i].mPEMapQuality = itpe->mPEMapQuality;
                itpe->mMerged = true;
            }
        }
    }
    util::loginfo("Finish augmenting SR supported SV candidates with DPs");
    // Append DP supported only SV 
    std::copy_if(pe.begin(), pe.end(), std::back_inserter(svs), [&](const SVRecord& sv){return !sv.mMerged;});
    pe.clear();
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
