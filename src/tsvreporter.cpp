#include "stats.h"

void Stats::reportTSV(const SVSet& svs, const GeneInfoList& gl){
    std::ofstream fw(mOpt->tsvOut);
    fw << "svType\tsvSize\tbpMark\tfuseGene\t";
    fw << "bp1Chr\tbp1Pos\tbp1Gene\t";
    fw << "bp2Chr\tbp2Pos\tbp2Gene\t";
    fw << "srCount\tdpCount\tsrRescued\tdpRescued\t";
    fw << "srRefCount\tdpRefCount\tAF\tinsBp\tinsSeq\t";
    fw << "bp1Trs\tbp2Trs\tsvSeq\tsvID\n";
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
        fw << gl[i].mFuseGene << "\t";
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
        // svID
        fw << svs[i].mID << "\n";
    }
    fw.close();
}
