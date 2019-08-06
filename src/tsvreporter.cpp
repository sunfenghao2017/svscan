#include "stats.h"

void Stats::reportTSV(const SVSet& svs, const GeneInfoList& gl){
    std::ofstream fw(mOpt->tsvOut);
    fw << "SVType\tCatType\t";
    fw << "Chr1\tBreakpoint1\tGene1\t";
    fw << "Chr2\tBreakpoing2\tGene2\t";
    fw << "SRSupports\tDPSupports\t";
    fw << "RefCount\tAltCount\tFreq\t";
    fw << "Strand1\tTranscripts1\t";
    fw << "Strand2\tTranscripts2\tSVID\n";
    for(uint32_t i = 0; i < gl.size(); ++i){
        fw << svutil::addID(gl[i].mSVT) << "\t" << svutil::addOrientation(gl[i].mSVT) << "\t";
        fw << gl[i].mChr1 << "\t" << gl[i].mPos1 << "\t" << gl[i].mGene1 << "\t";
        fw << gl[i].mChr2 << "\t" << gl[i].mPos2 << "\t" << gl[i].mGene2 << "\t";
        fw << mJctCnts[i].mAltQual.size() << "\t" << mSpnCnts[i].mAltQual.size() << "\t";
        if(svs[i].mPrecise){
            fw << mJctCnts[i].mRefQual.size() << "\t" << mJctCnts[i].mAltQual.size() << "\t";
            fw << (double)mJctCnts[i].mAltQual.size()/(double)(mJctCnts[i].mRefQual.size() + mJctCnts[i].mAltQual.size()) << "\t";
        }else{
            fw << mSpnCnts[i].mRefQual.size() << "\t" << mSpnCnts[i].mAltQual.size() << "\t";
            fw << (double)mSpnCnts[i].mAltQual.size()/(double)(mSpnCnts[i].mRefQual.size() + mSpnCnts[i].mAltQual.size()) << "\t";
        }
        fw << util::join(gl[i].mStrand1, ",") << "\t" << util::join(gl[i].mTrans1, ",") << "\t";
        fw << util::join(gl[i].mStrand2, ",") << "\t" << util::join(gl[i].mTrans2, ",") << "\t";
        fw << svs[i].mID << "\n";
    }
    fw.close();
}
