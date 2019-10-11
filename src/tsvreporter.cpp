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
    fw << "bp1Trs\tbp2Trs\tsvSeq\tseqBp\tID\tsvtInt\tfsMask";//[25,31]
    if(mOpt->rnamode){
        fw << "\tts1Name\tts1Pos\tts2Name\tts2Pos\n"; //[32,35]
    }else{
        fw << "\n";
    }
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
        if(mOpt->rnamode){
            // bp1Chr bp1Pos bp1Gene
            fw << gl[i].mChr1 << "\t" << gl[i].mPos1 << "\t" << gl[i].mGene1 << "\t";
            // bp2Chr bp2Pos bp2Gene
            fw << gl[i].mChr2 << "\t" << gl[i].mPos2 << "\t" << gl[i].mGene2 << "\t";
        }else{
            // bp1Chr bp1Pos bp1Gene
            fw << svs[i].mNameChr1 << "\t" << svs[i].mSVStart << "\t" << gl[i].mGene1 << "\t";
            // bp2Chr bp2Pos bp2Gene
            fw << svs[i].mNameChr2 << "\t" << svs[i].mSVEnd << "\t" << gl[i].mGene2 << "\t";
        }
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
        if(gl[i].mTrans1.empty()) fw << ".\t";
        else fw << util::join(gl[i].mTrans1, ";") << "\t";
        if(gl[i].mTrans2.empty()) fw << ".\t";
        else fw << util::join(gl[i].mTrans2, ";") << "\t";
        // svSeq seqBp
        if(svs[i].mSVT == 4) fw << svs[i].mInsSeq << "\t-\t";
        else if(svs[i].mPrecise) fw << svs[i].mConsensus << "\t" << svs[i].mGapCoord[0] << "\t";
        else fw << "-\t-\t";
        // svID svtInt fsMask
        fw << svs[i].mID << "\t" << svs[i].mSVT << "\t" << gl[i].mFuseGene.status;
        if(mOpt->rnamode){
            // ts1Name ts1Pos ts2Name ts2Pos
            fw << "\t" << svs[i].mNameChr1 << "\t" << svs[i].mSVStart << "\t" << svs[i].mNameChr2 << "\t" << svs[i].mSVEnd << "\n";
        }else{
            fw << "\n";
        }
    }
    fw.close();
}

void Stats::maskFuseRec(const SVSet& svs, GeneInfoList& gl){
    mOpt->fuseOpt->init();
    // mask hot gene status
    std::map<std::string, std::set<std::string>> fpairs;
    for(uint32_t i = 0; i < gl.size(); ++i){
        if((gl[i].mFuseGene.status & FUSION_FALLGENE)){
            if(mOpt->fuseOpt->hasWhiteGene(gl[i].mFuseGene.hgene, gl[i].mFuseGene.tgene)){
                gl[i].mFuseGene.status |= FUSION_FHOTGENE;
                if(mOpt->fuseOpt->matchHotDirec(gl[i].mFuseGene.hgene, gl[i].mFuseGene.tgene)){
                    gl[i].mFuseGene.status |= FUSION_FCOMMONHOTDIRECT;
                }
            }
            if(gl[i].mFuseGene.hgene != gl[i].mFuseGene.tgene){
                fpairs[gl[i].mFuseGene.hgene].insert(gl[i].mFuseGene.tgene);
            }else{
                gl[i].mFuseGene.status |= FUSION_FINSAMEGENE;
            }
        }
    }
    // mask fusion mirror pair
    for(uint32_t i = 0; i < gl.size(); ++i){
        if((!(gl[i].mFuseGene.status & FUSION_FALLGENE)) || (gl[i].mFuseGene.status & FUSION_FINSAMEGENE)) continue;
        std::string hg = gl[i].mFuseGene.hgene;
        std::string tg = gl[i].mFuseGene.tgene;
        auto titer = fpairs.find(tg);
        if(titer != fpairs.end() && titer->second.find(hg) != titer->second.end()){
            gl[i].mFuseGene.status |= FUSION_FMIRROR;
        }
    }
    // mask other status of fusions
    for(uint32_t i = 0; i < gl.size(); ++i){
        if(!(gl[i].mFuseGene.status & FUSION_FALLGENE)) continue;
        if(mOpt->fuseOpt->hasBlackGene(gl[i].mFuseGene.hgene, gl[i].mFuseGene.tgene)){
            gl[i].mFuseGene.status |= FUSION_FBLACKPAIR;
        }
        if(mOpt->fuseOpt->inBlackList(gl[i].mFuseGene.hgene, gl[i].mFuseGene.tgene)){
            gl[i].mFuseGene.status |= FUSION_FBLACKGENE;
        }
        if(!mOpt->fuseOpt->validSV(svs[i].mSVT, svs[i].mNameChr1, svs[i].mNameChr2, svs[i].mSVStart, svs[i].mSVEnd)){
            gl[i].mFuseGene.status |= FUSION_FFBG;
        }
        if(mOpt->fuseOpt->inWhiteList(gl[i].mFuseGene.hgene, gl[i].mFuseGene.tgene)){
            gl[i].mFuseGene.status |= FUSION_FINDB;
        }
        if(mOpt->fuseOpt->hasWhiteGene(gl[i].mFuseGene.hgene, gl[i].mFuseGene.tgene)){
            gl[i].mFuseGene.status |= FUSION_FHOTGENE;
        }
        if(gl[i].mFuseGene.status & FUSION_FINSAMEGENE){
            if(svutil::trsUnitIsNear(gl[i].mTrans1, gl[i].mTrans2, 1)){
                gl[i].mFuseGene.status |= FUSION_FTOOSMALLSIZE;
            }
        }
        if(svs[i].mSVT == 2 && mOpt->rnamode){
            gl[i].mFuseGene.status |= FUSION_FTOOSMALLSIZE;
        }
        if(svs[i].mPrecise){
            gl[i].mFuseGene.status |= FUSION_FPRECISE;
            if(svutil::simpleSeq(svs[i].mConsensus.substr(0, svs[i].mGapCoord[0])) ||
               svutil::simpleSeq(svs[i].mConsensus.substr(svs[i].mGapCoord[1]))){
                gl[i].mFuseGene.status |= FUSION_FLOWCOMPLEX;
            }
        }
        float af = 0.0;
        int32_t srv = mJctCnts[i].mAltQual.size();
        int32_t srr = mJctCnts[i].mRefQual.size();
        int32_t dpv = mSpnCnts[i].mAltQual.size();
        int32_t dpr = mSpnCnts[i].mRefQual.size();
        if(svs[i].mPrecise) af = (double)(srv)/(double)(srv + srr);
        else af = (double)(dpv)/(double)(dpv + dpr);
        if(gl[i].mFuseGene.status & FUSION_FINDB){// fusion in whitelist
            if((srv < mOpt->fuseOpt->mWhiteFilter.mMinSupport) && (dpv < mOpt->fuseOpt->mWhiteFilter.mMinSupport)){
                gl[i].mFuseGene.status |= FUSION_FLOWSUPPORT;
            }
            if(af < mOpt->fuseOpt->mWhiteFilter.mMinVAF){
                gl[i].mFuseGene.status |= FUSION_FLOWAF;
            }
            if((srv + srr) < mOpt->fuseOpt->mWhiteFilter.mMinDepth && ((dpr + dpv) < mOpt->fuseOpt->mWhiteFilter.mMinDepth)){
                gl[i].mFuseGene.status |= FUSION_FLOWDEPTH;
            }
            if((svs[i].mSVT != 2) && gl[i].mFuseGene.status & FUSION_FINSAMEGENE){
                if(svs[i].mSize < mOpt->fuseOpt->mWhiteFilter.mMinIntraGeneSVSize){
                    gl[i].mFuseGene.status |= FUSION_FTOOSMALLSIZE;
                }
            }
        }else if(gl[i].mFuseGene.status & FUSION_FHOTGENE){// fusion not in whitelist
            if((srv < mOpt->fuseOpt->mUsualFilter.mMinSupport) && (dpv < mOpt->fuseOpt->mUsualFilter.mMinSupport)){
                gl[i].mFuseGene.status |= FUSION_FLOWSUPPORT;
            }
            if(af < mOpt->fuseOpt->mUsualFilter.mMinVAF){
                gl[i].mFuseGene.status |= FUSION_FLOWAF;
            }
            if((srv + srr) < mOpt->fuseOpt->mUsualFilter.mMinDepth && ((dpr + dpv) < mOpt->fuseOpt->mUsualFilter.mMinDepth)){
                gl[i].mFuseGene.status |= FUSION_FLOWDEPTH;
            }
            if((svs[i].mSVT != 2) && gl[i].mFuseGene.status & FUSION_FINSAMEGENE){
                if(svs[i].mSize < mOpt->fuseOpt->mUsualFilter.mMinIntraGeneSVSize){
                    gl[i].mFuseGene.status |= FUSION_FTOOSMALLSIZE;
                }
            }
        }
    }
    // mask primary and supplementary
    uint32_t REPORT_REQUEST = (FUSION_FINDB | FUSION_FHOTGENE);
    uint32_t ALL_DROP_MASK = (FUSION_FBLACKGENE | FUSION_FBLACKPAIR | FUSION_FFBG);
    uint32_t PRIMARY_DROP_MASK = (FUSION_FLOWAF | FUSION_FLOWSUPPORT | FUSION_FLOWDEPTH | FUSION_FLOWCOMPLEX | FUSION_FTOOSMALLSIZE);
    for(uint32_t i = 0; i < gl.size(); ++i){
        if(!(gl[i].mFuseGene.status & FUSION_FALLGENE)) continue;
        if(gl[i].mFuseGene.status & ALL_DROP_MASK) continue;
        if(gl[i].mFuseGene.status & REPORT_REQUEST){
            if((gl[i].mFuseGene.status & FUSION_FNORMALCATDIRECT) && 
               (gl[i].mFuseGene.status & FUSION_FCOMMONHOTDIRECT)){
                if(!(gl[i].mFuseGene.status & PRIMARY_DROP_MASK)){
                    gl[i].mFuseGene.status |= FUSION_FPRIMARYR;
                }
            }else{
                if(!(gl[i].mFuseGene.status & PRIMARY_DROP_MASK)){
                    gl[i].mFuseGene.status |= FUSION_FSUPPLEMENTARY;
                }
            }
        }
    }
}

void Stats::reportFusionTSV(const SVSet& svs, GeneInfoList& gl){
    // output valid fusions
    std::string header = "FusionGene\tFusionPattern\tFusionReads\tTotalReads\tFusionRate\t"; //[0-4]
    header.append("Gene1\tChr1\tJunctionPosition1\tStrand1\tTranscript1\t");//[5-9]
    header.append("Gene2\tChr2\tJunctionPosition2\tStrand2\tTranscript2\t");//[10-14]
    header.append("FusionSequence\tfseqBp\tinDB\tsvType\tsvSize\t"); //[15-19]
    header.append("srCount\tdpCount\tsrRescued\tdpRescued\tsrRefCount\tdpRefCount\t"); //[20-25]
    header.append("insBp\tinsSeq\tsvID\tsvtInt\tfsMask"); //[26-29]
    if(mOpt->rnamode){
        header.append("\tts1Name\tts1Pos\tts2Name\tts2Pos\n"); //[30-33]
    }else{
        header.append("\n");
    }
    std::ofstream fw(mOpt->fuseOpt->mOutFile);
    std::ofstream fs(mOpt->fuseOpt->mSupFile);
    fw << header;
    fs << header;
    for(uint32_t i = 0; i < gl.size(); ++i){
        if(gl[i].mFuseGene.status & FUSION_FPRIMARYR){
            fw << toFuseRec(svs, gl, i);
        }
        if(gl[i].mFuseGene.status & FUSION_FSUPPLEMENTARY){
            fs << toFuseRec(svs, gl, i);
        }
    }
    fw.close();
    fs.close();
}

std::string Stats::toFuseRec(const SVSet& svs, const GeneInfoList& gl, int32_t i){
    std::stringstream oss;
    float af = 0.0;
    int32_t srv = mJctCnts[i].mAltQual.size();
    int32_t srr = mJctCnts[i].mRefQual.size();
    int32_t dpv = mSpnCnts[i].mAltQual.size();
    int32_t dpr = mSpnCnts[i].mRefQual.size();
    if(svs[i].mPrecise) af = (double)(srv)/(double)(srv + srr);
    else af = (double)(dpv)/(double)(dpv + dpr);
    oss << gl[i].mFuseGene.hgene << "->" << gl[i].mFuseGene.tgene << "\t"; // FusionGene
    if(gl[i].mGene1 == gl[i].mFuseGene.hgene){// FusionPattern
        oss << gl[i].mStrand1 << gl[i].mStrand2 << "\t";
    }else{
        oss << gl[i].mStrand2 << gl[i].mStrand1 << "\t";
    }
    if(svs[i].mPrecise){// FusionReads TotalReads FusionRate
        oss << srv << "\t";
        oss << (srr + srv) << "\t"; 
        oss << af << "\t";
    }else{
        oss << dpv << "\t";
        oss << (dpr + dpv) << "\t";
        oss << af << "\t";
    }
    if(gl[i].mGene1 == gl[i].mFuseGene.hgene){
        // Gene1 Chr1 JunctionPosition1 Strand1 Transcript1
        oss << gl[i].mGene1 << "\t";
        if(mOpt->rnamode){
            oss << gl[i].mChr1 << "\t" << gl[i].mPos1 << "\t";
        }else{
            oss << svs[i].mNameChr1 << "\t" << svs[i].mSVStart << "\t";
        }
        oss <<  gl[i].mStrand1 << "\t"  << util::join(gl[i].mTrans1, ";") << "\t";
        // Gene2 Chr2 JunctionPosition2 Strand2 Transcript2
        oss << gl[i].mGene2 << "\t";
        if(mOpt->rnamode){
            oss << gl[i].mChr2 << "\t" << gl[i].mPos2 << "\t";
        }else{
            oss << svs[i].mNameChr2 << "\t" << svs[i].mSVStart << "\t";
        }
        oss <<  gl[i].mStrand2 << "\t"  << util::join(gl[i].mTrans2, ";") << "\t";
    }else{
        // Gene2 Chr2 JunctionPosition2 Strand2 Transcript2
        if(mOpt->rnamode){
            oss << gl[i].mChr2 << "\t" << gl[i].mPos2 << "\t";
        }else{
            oss << svs[i].mNameChr2 << "\t" << svs[i].mSVStart << "\t";
        }
        oss <<  gl[i].mStrand2 << "\t"  << util::join(gl[i].mTrans2, ";") << "\t";
        // Gene1 Chr1 JunctionPosition1 Strand1 Transcript1
        if(mOpt->rnamode){
            oss << gl[i].mChr1 << "\t" << gl[i].mPos1 << "\t";
        }else{
            oss << svs[i].mNameChr1 << "\t" << svs[i].mSVStart << "\t";
        }
        oss <<  gl[i].mStrand1 << "\t"  << util::join(gl[i].mTrans1, ";") << "\t";
    }
    // FusinSequence
    if(svs[i].mSVT == 4) oss << svs[i].mInsSeq << "\t-\t";
    else if(svs[i].mPrecise) oss << svs[i].mConsensus << "\t" << svs[i].mGapCoord[0] << "\t";
    else oss << "-\t-\t";
    if(gl[i].mFuseGene.status & FUSION_FINDB) oss << "Y\t"; // inDB
    else oss << "N\t";
    oss << svutil::addID(svs[i].mSVT) << "\t";  // svType
    if(svs[i].mSVT >= 5) oss << "-\t";
    else oss << svs[i].mSize << "\t";           // svSize
    oss << svs[i].mSRSupport << "\t";           // srCount
    oss << svs[i].mPESupport << "\t";           // dpCount
    oss << mJctCnts[i].mAltQual.size() << "\t"; // srRescued
    oss << mSpnCnts[i].mAltQual.size() << "\t"; // dpRescued
    oss << mJctCnts[i].mRefQual.size() << "\t"; // srRefCount
    oss <<mSpnCnts[i].mRefQual.size() << "\t";  // dpRefCount
    oss << svs[i].mBpInsSeq.length() << "\t";   // insBp
    if(svs[i].mBpInsSeq.empty()) oss << "-\t";
    else oss << svs[i].mBpInsSeq << "\t";       // insSeq
    oss << svs[i].mID << "\t";                  // svID
    oss << svs[i].mSVT << "\t";                 // svtInt
    oss << gl[i].mFuseGene.status;              // fsMask
    if(mOpt->rnamode){
        oss << "\t";
        oss << svs[i].mNameChr1 << "\t";            // ts1Name
        oss << svs[i].mSVStart << "\t";         // ts2Pos
        oss << svs[i].mNameChr2 << "\t";            // ts2Name
        oss << svs[i].mSVEnd << "\n";           // ts2Pos
    }else{
        oss << "\n";
    }
    return oss.str();
}
