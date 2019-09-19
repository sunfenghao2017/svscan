#include "stats.h"
#include "fusionopt.h"

void Stats::reportSVTSV(const SVSet& svs, const GeneInfoList& gl){
    std::ofstream fw(mOpt->tsvOut);
    fw << "svType\tsvSize\tbpMark\t";//[0, 2]
    fw << "fuseGene\t"; //[3];
    fw << "hGene\thTrsEnd\thTrsStrand\t";//[4,6]
    fw << "tGene\ttTrsEnd\ttTrsStrand\t";//[7,9]
    fw << "bp1Chr\tbp1Pos\tbp1Gene\t";//[10,12]
    fw << "bp2Chr\tbp2Pos\tbp2Gene\t";//[13,15]
    fw << "srCount\tdpCount\tsrRescued\tdpRescued\t";//[16,19]
    fw << "srRefCount\tdpRefCount\tAF\tinsBp\tinsSeq\t";//[20,24]
    fw << "bp1Trs\tbp2Trs\tsvSeq\tsvID\tsvtInt\n";//[25,29]
    for(uint32_t i = 0; i < gl.size(); ++i){
        // skip false positive insertions
        if(svs[i].mSVT == 4 && mJctCnts[i].mFPIns > mSpnCnts[i].mAltQual.size() * mOpt->filterOpt->mMaxFPIns) continue;
        // svType
        fw << svutil::addID(svs[i].mSVT) << "\t";
        // svSize
        if(svs[i].mSVT >= 5) fw << "-" << "\t";
        else if(svs[i].mSVT == 4) fw << svs[i].mInsSeq.size() << "\t";
        else fw << svs[i].mSize << "\t";
        // bpMark
        fw << svutil::getBpMark(svs[i].mSVT) << "\t";
        // fuseGene
        fw << gl[i].mFuseGene;
        // bp1Chr bp1Pos bp1Gene
        fw << svs[i].mNameChr1 << "\t" << svs[i].mSVStart << "\t" << gl[i].mGene1 << "\t";
        // bp2Chr bp2Pos bp2Gene
        fw << svs[i].mNameChr2 << "\t" << svs[i].mSVEnd << "\t" << gl[i].mGene2 << "\t";
        // srCount dpCount srRescued dpRescued
        fw << svs[i].mSRSupport << "\t" << svs[i].mPESupport << "\t" << mJctCnts[i].mAltQual.size() << "\t" << mSpnCnts[i].mAltQual.size() << "\t";
        // srRefCount dpRefCount
        fw << mJctCnts[i].mRefQual.size() << "\t" <<mSpnCnts[i].mRefQual.size() << "\t";
        // AF
        if(svs[i].mPrecise) fw << (double)(mJctCnts[i].mAltQual.size())/(double)(mJctCnts[i].mRefQual.size() + mJctCnts[i].mAltQual.size()) << "\t";
        else fw << (double)(mSpnCnts[i].mAltQual.size())/(double)(mSpnCnts[i].mRefQual.size() + mSpnCnts[i].mAltQual.size()) << "\t";
        // insBp insSeq
        fw << svs[i].mBpInsSeq.length() << "\t" << (svs[i].mBpInsSeq.length() == 0 ? "-" : svs[i].mBpInsSeq) << "\t"; 
        // bp1Trs bp2Trs svID
        fw << util::join(gl[i].mTrans1, ",") << "\t";
        fw << util::join(gl[i].mTrans2, ",") << "\t";
        // svSeq
        if(svs[i].mSVT == 4) fw << svs[i].mInsSeq << "\t";
        else if(svs[i].mPrecise) fw << svs[i].mConsensus << "\t";
        else fw << "-" << "\t";
        // svID svtInt
        fw << svs[i].mID << "\t" << svs[i].mSVT << "\n";

    }
    fw.close();
}

void Stats::reportFusionTSV(const SVSet& svs, const GeneInfoList& gl){
    mOpt->fuseOpt->init();
    std::ofstream fw(mOpt->fuseOpt->mOutFile);
    fw << "FusionGene\tFusionPattern\tFusionReads\tTotalReads\tFusionRate\t"; //[0-4]
    fw << "Gene1\tChr1\tJunctionPosition1\tStrand1\tTranscript1\t";//[5-9]
    fw << "Gene2\tChr2\tJunctionPosition2\tStrand2\tTranscript2\t";//[10-14]
    fw << "FusionSequence\tinDB\tsvType\tsvSize\t"; //[15-17]
    fw << "srCount\tdpCount\tsrRescued\tdpRescued\tsrRefCount\tdpRefCount\t"; //[18-23]
    fw << "insBp\tinsSeq\tsvID\tsvtInt\n"; //[24-27]
    for(uint32_t i = 0; i < gl.size(); ++i){
        // keep only (hgene+5'->tgene+3') fusion
        if(!gl[i].mFuseGene.valid) continue;
        // skip fusion which does not contain any gene in whitelist
        if(!mOpt->fuseOpt->hasWhiteGene(gl[i].mFuseGene.hgene, gl[i].mFuseGene.tgene)) continue;
        // skip fusion in blacklist
        if(mOpt->fuseOpt->inBlackList(gl[i].mFuseGene.hgene, gl[i].mFuseGene.tgene)) continue;
        bool inWhitelist = mOpt->fuseOpt->inWhiteList(gl[i].mFuseGene.hgene, gl[i].mFuseGene.tgene);
        bool keep = false;
        float af = 0.0;
        int32_t srv = mJctCnts[i].mAltQual.size();
        int32_t srr = mJctCnts[i].mRefQual.size();
        int32_t dpv = mSpnCnts[i].mAltQual.size();
        int32_t dpr = mSpnCnts[i].mRefQual.size();
        if(svs[i].mPrecise) af = (double)(srv)/(double)(srv + srr);
        else af = (double)(dpv)/(double)(dpv + dpr);
        if(inWhitelist){// fusion in whitelist
           if((srv >= mOpt->fuseOpt->mWhiteFilter.mMinSupport || dpv >= mOpt->fuseOpt->mWhiteFilter.mMinSupport) &&
              (af >= mOpt->fuseOpt->mWhiteFilter.mMinVAF)) keep = true;
           if((srv + srr) < mOpt->fuseOpt->mWhiteFilter.mMinDepth && ((dpr + dpv) < mOpt->fuseOpt->mWhiteFilter.mMinDepth)){
               keep = false;
           }
        }else{// fusion not in whitelist
            if((srv >= mOpt->fuseOpt->mUsualFilter.mMinSupport || dpv >= mOpt->fuseOpt->mUsualFilter.mMinSupport) &&
               (af >= mOpt->fuseOpt->mUsualFilter.mMinVAF)) keep = true;
            if((srv + srr) < mOpt->fuseOpt->mUsualFilter.mMinDepth && ((dpr + dpv) < mOpt->fuseOpt->mUsualFilter.mMinDepth)){
               keep = false;
           }
        }
        if(!keep) continue;
        // skip fusion in background
        if(!mOpt->fuseOpt->validSV(svs[i].mSVT, svs[i].mNameChr1, svs[i].mNameChr2, svs[i].mSVStart, svs[i].mSVEnd)) continue;
        fw << gl[i].mFuseGene.hgene << "->" << gl[i].mFuseGene.tgene << "\t"; // FusionGene
        if(gl[i].mGene1 == gl[i].mFuseGene.hgene){// FusionPattern
            fw << gl[i].mStrand1 << gl[i].mStrand2 << "\t";
        }else{
            fw << gl[i].mStrand2 << gl[i].mStrand1 << "\t";
        }
        if(svs[i].mPrecise){// FusionReads TotalReads FusionRate
            fw << srv << "\t";
            fw << (srr + srv) << "\t"; 
            fw << af << "\t";
        }else{
            fw << dpv << "\t";
            fw << (dpr + dpv) << "\t";
            fw << af << "\t";
        }
        if(gl[i].mGene1 == gl[i].mFuseGene.hgene){
            // Gene1 Chr1 JunctionPosition1 Strand1 Transcript1
            fw << gl[i].mGene1 << "\t" << svs[i].mNameChr1 << "\t" << svs[i].mSVStart << "\t" <<  gl[i].mStrand1 << "\t"  << util::join(gl[i].mTrans1, ",") << "\t";
            // Gene2 Chr2 JunctionPosition2 Strand2 Transcript2
            fw << gl[i].mGene2 << "\t" << svs[i].mNameChr2 << "\t" << svs[i].mSVEnd << "\t" <<  gl[i].mStrand2 << "\t"  << util::join(gl[i].mTrans2, ",") << "\t";
        }else{
            // Gene2 Chr2 JunctionPosition2 Strand2 Transcript2
            fw << gl[i].mGene2 << "\t" << svs[i].mNameChr2 << "\t" << svs[i].mSVEnd << "\t" <<  gl[i].mStrand2 << "\t"  << util::join(gl[i].mTrans2, ",") << "\t";
            // Gene1 Chr1 JunctionPosition1 Strand1 Transcript1
            fw << gl[i].mGene1 << "\t" << svs[i].mNameChr1 << "\t" << svs[i].mSVStart << "\t" <<  gl[i].mStrand1 << "\t"  << util::join(gl[i].mTrans1, ",") << "\t";
        }
        // FusinSequence
        if(svs[i].mConsensus.empty()) fw << "-\t";
        else fw << svs[i].mConsensus << "\t";
        if(inWhitelist) fw << "Y\t";               // inDB
        else fw << "N\t";
        fw << svutil::addID(svs[i].mSVT) << "\t";  // svType
        if(svs[i].mSVT >= 5) fw << "-\t";
        else fw << svs[i].mSize << "\t";           // svSize
        fw << svs[i].mSRSupport << "\t";           // srCount
        fw << svs[i].mPESupport << "\t";           // dpCount
        fw << mJctCnts[i].mAltQual.size() << "\t"; // srRescued
        fw << mSpnCnts[i].mAltQual.size() << "\t"; // dpRescued
        fw << mJctCnts[i].mRefQual.size() << "\t"; // srRefCount
        fw <<mSpnCnts[i].mRefQual.size() << "\t";  // dpRefCount
        fw << svs[i].mBpInsSeq.length() << "\t";   // insBp
        if(svs[i].mBpInsSeq.empty()) fw << "-\t";
        else fw << svs[i].mBpInsSeq << "\t";       // insSeq
        fw << svs[i].mID << "\t";                  // svID
        fw << svs[i].mSVT << "\n";                 // svtInt
    }
    fw.close();
}
