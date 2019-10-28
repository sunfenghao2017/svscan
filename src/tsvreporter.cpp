#include "stats.h"
#include "fusionopt.h"

void Stats::reportSVTSV(SVSet& svs, GeneInfoList& gl){
    std::ofstream fw(mOpt->tsvOut);
    fw << "svType\tsvSize\tbpMark\t";//[0,2]
    fw << "bp1Chr\tbp1Pos\tbp2Chr\tbp2Pos\t"; // [3,6]
    fw << "srCount\tdpCount\t";// [7,8]
    fw << "srRescued\tdpRescued\t"; // [9,10]
    fw << "srRefCount\tdpRefCount\tAF\t"; // [11,13]
    fw << "insBp\tinsSeq\tsvSeq\tseqBp\t";// [14,17]
    fw << "ID\tsvtInt\t"; // [18,19]
    fw << "bp1Gene\tbp2Gene\tfuseGene\tfsMask"; // [20,23]
    if(mOpt->rnamode){
        fw << "\tts1Name\tts1Pos\tts2Name\tts2Pos\n"; //[24,27]
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
        // bp1Chr bp1Pos bp2Chr bp2Pos
        fw << gl[i].mChr1 << "\t" << gl[i].mPos1 << "\t" << gl[i].mChr2 << "\t" << gl[i].mPos2 << "\t";
        // srCount dpCount srRescued dpRescued
        fw << svs[i].mSRSupport << "\t" << svs[i].mPESupport << "\t" << mJctCnts[i].mAltQual.size() << "\t" << mSpnCnts[i].mAltQual.size() << "\t";
        // srRefCount dpRefCount
        fw << mJctCnts[i].mRefQual.size() << "\t" <<mSpnCnts[i].mRefQual.size() << "\t";
        // AF
        if(svs[i].mPrecise) fw << (double)(mJctCnts[i].mAltQual.size())/(double)(mJctCnts[i].mRefQual.size() + mJctCnts[i].mAltQual.size()) << "\t";
        else fw << (double)(mSpnCnts[i].mAltQual.size())/(double)(mSpnCnts[i].mRefQual.size() + mSpnCnts[i].mAltQual.size()) << "\t";
        // insBp insSeq
        fw << svs[i].mBpInsSeq.length() << "\t" << (svs[i].mBpInsSeq.length() == 0 ? "-" : svs[i].mBpInsSeq) << "\t"; 
        // svSeq seqBp
        if(svs[i].mSVT == 4) fw << svs[i].mInsSeq << "\t-\t";
        else if(svs[i].mPrecise) fw << svs[i].mConsensus << "\t" << svs[i].mGapCoord[0] << "\t";
        else fw << "-\t-\t";
        // svID svtInt
        fw << svs[i].mID << "\t" << svs[i].mSVT << "\t";
        // bp1Gene bp2Gene svID
        if(gl[i].mGene1.empty()) fw << "-\t";
        else fw << gl[i].getTrs1() << "\t";
        if(gl[i].mGene2.empty()) fw << "-\t";
        else fw << gl[i].getTrs2() << "\t";
        // fuseGene fsMask
        fw << gl[i].getFuseGene() << "\t" << gl[i].getFsMask();
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
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            if(mOpt->rnamode) gl[i].mFuseGene[j].status |= FUSION_FCALLFROMRNASEQ; // mask rna/dna calling
            if((gl[i].mFuseGene[j].status & FUSION_FALLGENE)){
                if(mOpt->fuseOpt->hasWhiteGene(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                    gl[i].mFuseGene[j].status |= FUSION_FHOTGENE;
                    if(mOpt->fuseOpt->matchHotDirec(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                        gl[i].mFuseGene[j].status |= FUSION_FCOMMONHOTDIRECT;
                    }
                }
                if(gl[i].mFuseGene[j].hgene != gl[i].mFuseGene[j].tgene){
                    fpairs[gl[i].mFuseGene[j].hgene].insert(gl[i].mFuseGene[j].tgene);
                }else{
                    gl[i].mFuseGene[j].status |= FUSION_FINSAMEGENE;
                }
            }
        }
    }
    // mask fusion mirror pair
    for(uint32_t i = 0; i < gl.size(); ++i){
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            if((!(gl[i].mFuseGene[j].status & FUSION_FALLGENE)) || (gl[i].mFuseGene[j].status & FUSION_FINSAMEGENE)) continue;
            std::string hg = gl[i].mFuseGene[j].hgene;
            std::string tg = gl[i].mFuseGene[j].tgene;
            auto titer = fpairs.find(tg);
            if(titer != fpairs.end() && titer->second.find(hg) != titer->second.end()){
                gl[i].mFuseGene[j].status |= FUSION_FMIRROR;
            }
        }
    }
    // mask other status of fusions
    for(uint32_t i = 0; i < gl.size(); ++i){
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            if(!(gl[i].mFuseGene[j].status & FUSION_FALLGENE)) continue;
            if(mOpt->fuseOpt->hasBlackGene(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                gl[i].mFuseGene[j].status |= FUSION_FBLACKGENE;
            }
            if(mOpt->fuseOpt->inBlackList(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                gl[i].mFuseGene[j].status |= FUSION_FBLACKPAIR;
            }
            if(!mOpt->fuseOpt->validSV(svs[i].mSVT, svs[i].mNameChr1, svs[i].mNameChr2, svs[i].mSVStart, svs[i].mSVEnd)){
                gl[i].mFuseGene[j].status |= FUSION_FFBG;
            }
            if(mOpt->fuseOpt->inWhiteList(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                gl[i].mFuseGene[j].status |= FUSION_FINDB;
            }
            if(mOpt->fuseOpt->inWhiteList(gl[i].mFuseGene[j].tgene, gl[i].mFuseGene[j].hgene)){
                gl[i].mFuseGene[j].status |= FUSION_FMIRRORINDB;
            }
            if(mOpt->fuseOpt->hasWhiteGene(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                gl[i].mFuseGene[j].status |= FUSION_FHOTGENE;
            }
            if(gl[i].mFuseGene[j].status & FUSION_FINSAMEGENE){
                if(svutil::trsUnitIsNear(gl[i].getTrs1(), gl[i].getTrs2(), 1)){
                    gl[i].mFuseGene[j].status |= FUSION_FTOOSMALLSIZE;
                }
            }
            if(svs[i].mPrecise){
                gl[i].mFuseGene[j].status |= FUSION_FPRECISE;
                if(svutil::simpleSeq(svs[i].mConsensus.substr(0, svs[i].mGapCoord[0])) ||
                   svutil::simpleSeq(svs[i].mConsensus.substr(svs[i].mGapCoord[1]))){
                    gl[i].mFuseGene[j].status |= FUSION_FLOWCOMPLEX;
                }
            }
            float af = 0.0;
            int32_t srv = mJctCnts[i].mAltQual.size();
            int32_t srr = mJctCnts[i].mRefQual.size();
            int32_t dpv = mSpnCnts[i].mAltQual.size();
            int32_t dpr = mSpnCnts[i].mRefQual.size();
            if(svs[i].mPrecise) af = (double)(srv)/(double)(srv + srr);
            else af = (double)(dpv)/(double)(dpv + dpr);
            if(gl[i].mFuseGene[j].status & (FUSION_FINDB | FUSION_FMIRRORINDB)){// fusion in public database
                if((srv < mOpt->fuseOpt->mWhiteFilter.mMinSupport) && (dpv < mOpt->fuseOpt->mWhiteFilter.mMinSupport)){
                    gl[i].mFuseGene[j].status |= FUSION_FLOWSUPPORT;
                }
                if(af < mOpt->fuseOpt->mWhiteFilter.mMinVAF){
                    gl[i].mFuseGene[j].status |= FUSION_FLOWAF;
                }
                if((srv + srr) < mOpt->fuseOpt->mWhiteFilter.mMinDepth && ((dpr + dpv) < mOpt->fuseOpt->mWhiteFilter.mMinDepth)){
                    gl[i].mFuseGene[j].status |= FUSION_FLOWDEPTH;
                }
                if((svs[i].mSVT != 2) && gl[i].mFuseGene[j].status & FUSION_FINSAMEGENE){
                    if(svs[i].mSize < mOpt->fuseOpt->mWhiteFilter.mMinIntraGeneSVSize){
                        gl[i].mFuseGene[j].status |= FUSION_FTOOSMALLSIZE;
                    }
                }
            }else if(gl[i].mFuseGene[j].status & FUSION_FHOTGENE){// fusion in whitelist
                if((srv < mOpt->fuseOpt->mUsualFilter.mMinSupport) && (dpv < mOpt->fuseOpt->mUsualFilter.mMinSupport)){
                    gl[i].mFuseGene[j].status |= FUSION_FLOWSUPPORT;
                }
                if(af < mOpt->fuseOpt->mUsualFilter.mMinVAF){
                    gl[i].mFuseGene[j].status |= FUSION_FLOWAF;
                }
                if((srv + srr) < mOpt->fuseOpt->mUsualFilter.mMinDepth && ((dpr + dpv) < mOpt->fuseOpt->mUsualFilter.mMinDepth)){
                    gl[i].mFuseGene[j].status |= FUSION_FLOWDEPTH;
                }
                if((svs[i].mSVT != 2) && gl[i].mFuseGene[j].status & FUSION_FINSAMEGENE){
                    if(svs[i].mSize < mOpt->fuseOpt->mUsualFilter.mMinIntraGeneSVSize){
                        gl[i].mFuseGene[j].status |= FUSION_FTOOSMALLSIZE;
                    }
                }
            }
        }
    }
    // drop bits mask of all fusion events, if an fusion match any bit in FUSION_DROP_MASK, it will not be reported
    TFUSION_FLAG FUSION_DROP_MASK = (FUSION_FBLACKGENE | FUSION_FBLACKPAIR  | FUSION_FFBG | FUSION_FLOWCOMPLEX |
                                     FUSION_FTOOSMALLSIZE | FUSION_FLOWAF | FUSION_FLOWDEPTH);
    if(mOpt->rnamode) FUSION_DROP_MASK |= FUSION_FINSAMEGENE;
    // primary keep bits mask, fusion reported as primary must match all the bits in PRIMARY_KEEP_MASK
    TFUSION_FLAG PRIMARY_KEEP_MASK = (FUSION_FNORMALCATDIRECT | FUSION_FCOMMONHOTDIRECT | FUSION_FINDB);
    // keep bits mask, an fusion to be reported must match all bits in FUSION_KEEP_MASK
    TFUSION_FLAG FUSION_KEEP_MASK = (FUSION_FALLGENE | FUSION_FHOTGENE);
    for(uint32_t i = 0; i < gl.size(); ++i){
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            if(gl[i].mFuseGene[j].status & FUSION_DROP_MASK) continue;
            if((gl[i].mFuseGene[j].status & FUSION_KEEP_MASK) != FUSION_KEEP_MASK) continue;
            if((gl[i].mFuseGene[j].status & PRIMARY_KEEP_MASK) == PRIMARY_KEEP_MASK){
                gl[i].mFuseGene[j].status |= FUSION_FPRIMARY;
            }else{
                gl[i].mFuseGene[j].status |= FUSION_FSUPPLEMENTARY;
            }
        }
    }
}

void Stats::reportFusionTSV(SVSet& svs, GeneInfoList& gl){
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
        bool reported = false;
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            if(gl[i].mFuseGene[j].status & FUSION_FPRIMARY){
                fw << toFuseRec(svs[i], gl[i], j);
                reported = true;
            }
        }
        if(!reported){
            for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
                if(gl[i].mFuseGene[j].status & FUSION_FSUPPLEMENTARY){
                    fs << toFuseRec(svs[i], gl[i], j);
                }
            }
        }
    }
    fw.close();
    fs.close();
}

std::string Stats::toFuseRec(SVRecord& svr, GeneInfo& gi, int32_t i){
    std::stringstream oss;
    float af = 0.0;
    int32_t srv = mJctCnts[svr.mID].mAltQual.size();
    int32_t srr = mJctCnts[svr.mID].mRefQual.size();
    int32_t dpv = mSpnCnts[svr.mID].mAltQual.size();
    int32_t dpr = mSpnCnts[svr.mID].mRefQual.size();
    if(svr.mPrecise) af = (double)(srv)/(double)(srv + srr);
    else af = (double)(dpv)/(double)(dpv + dpr);
    oss << gi.mFuseGene[i].hgene << "->" << gi.mFuseGene[i].tgene << "\t"; // FusionGene
    // FusionPattern
    if(gi.mFuseGene[i].hfrom1) oss << gi.mGene1[gi.mFuseGene[i].hidx].strand;
    else oss << gi.mGene2[gi.mFuseGene[i].hidx].strand;
    if(gi.mFuseGene[i].tfrom1) oss << gi.mGene1[gi.mFuseGene[i].tidx].strand << "\t";
    else oss << gi.mGene2[gi.mFuseGene[i].tidx].strand << "\t";
    if(svr.mPrecise){// FusionReads TotalReads FusionRate
        oss << srv << "\t";
        oss << (srr + srv) << "\t"; 
        oss << af << "\t";
    }else{
        oss << dpv << "\t";
        oss << (dpr + dpv) << "\t";
        oss << af << "\t";
    }
    if(gi.mGene1[i].gene == gi.mFuseGene[i].hgene){
        // Gene1 Chr1 JunctionPosition1 Strand1 Transcript1
        oss << gi.mGene1[gi.mFuseGene[i].hidx].gene << "\t";
        if(mOpt->rnamode){
            oss << gi.mChr1 << "\t" << gi.mPos1 << "\t";
        }else{
            oss << svr.mNameChr1 << "\t" << svr.mSVStart << "\t";
        }
        oss <<  gi.mGene1[gi.mFuseGene[i].hidx].strand << "\t";
        oss << gi.mGene1[gi.mFuseGene[i].hidx].name << "\t";
        // Gene2 Chr2 JunctionPosition2 Strand2 Transcript2
        oss << gi.mGene2[gi.mFuseGene[i].tidx].gene << "\t";
        if(mOpt->rnamode){
            oss << gi.mChr2 << "\t" << gi.mPos2 << "\t";
        }else{
            oss << svr.mNameChr2 << "\t" << svr.mSVEnd << "\t";
        }
        oss <<  gi.mGene2[gi.mFuseGene[i].tidx].strand << "\t";
        oss << gi.mGene2[gi.mFuseGene[i].tidx].name << "\t";
    }else{
        // Gene1 Chr1 JunctionPosition1 Strand1 Transcript1
        oss << gi.mGene2[gi.mFuseGene[i].hidx].gene << "\t";
        if(mOpt->rnamode){
            oss << gi.mChr2 << "\t" << gi.mPos2 << "\t";
        }else{
            oss << svr.mNameChr2 << "\t" << svr.mSVStart << "\t";
        }
        oss <<  gi.mGene2[gi.mFuseGene[i].hidx].strand << "\t";
        oss << gi.mGene2[gi.mFuseGene[i].hidx].name << "\t";
        // Gene2 Chr2 JunctionPosition2 Strand2 Transcript2
        oss << gi.mGene1[gi.mFuseGene[i].tidx].gene << "\t";
        if(mOpt->rnamode){
            oss << gi.mChr1 << "\t" << gi.mPos1 << "\t";
        }else{
            oss << svr.mNameChr1 << "\t" << svr.mSVEnd << "\t";
        }
        oss <<  gi.mGene1[gi.mFuseGene[i].tidx].strand << "\t";
        oss << gi.mGene1[gi.mFuseGene[i].tidx].name << "\t";
    }
    // FusinSequence
    if(svr.mSVT == 4) oss << svr.mInsSeq << "\t-\t";
    else if(svr.mPrecise) oss << svr.mConsensus << "\t" << svr.mGapCoord[0] << "\t";
    else oss << "-\t-\t";
    if(gi.mFuseGene[i].status & FUSION_FINDB) oss << "Y\t"; // inDB
    else oss << "N\t";
    oss << svutil::addID(svr.mSVT) << "\t";                 // svType
    if(svr.mSVT >= 5) oss << "-\t";
    else oss << svr.mSize << "\t";                          // svSize
    oss << svr.mSRSupport << "\t";                          // srCount
    oss << svr.mPESupport << "\t";                          // dpCount
    oss << mJctCnts[svr.mID].mAltQual.size() << "\t";       // srRescued
    oss << mSpnCnts[svr.mID].mAltQual.size() << "\t";       // dpRescued
    oss << mJctCnts[svr.mID].mRefQual.size() << "\t";       // srRefCount
    oss <<mSpnCnts[svr.mID].mRefQual.size() << "\t";        // dpRefCount
    oss << svr.mBpInsSeq.length() << "\t";                  // insBp
    if(svr.mBpInsSeq.empty()) oss << "-\t";
    else oss << svr.mBpInsSeq << "\t";                      // insSeq
    oss << svr.mID << "\t";                                 // svID
    oss << svr.mSVT << "\t";                                // svtInt
    oss << gi.mFuseGene[i].status;                          // fsMask
    if(mOpt->rnamode){
        if(gi.mGene1[i].gene == gi.mFuseGene[i].hgene){
            oss << "\t";
            oss << svr.mNameChr1 << "\t";                   // ts1Name
            oss << svr.mSVStart << "\t";                    // ts1Pos
            oss << svr.mNameChr2 << "\t";                   // ts2Name
            oss << svr.mSVEnd << "\n";                      // ts2Pos
        }else{
            oss << "\t";
            oss << svr.mNameChr2 << "\t";                   // ts1Name
            oss << svr.mSVEnd << "\t";                      // ts1Pos
            oss << svr.mNameChr1 << "\t";                   // ts2Name
            oss << svr.mSVStart << "\n";                    // ts2Pos
        }
    }else{
        oss << "\n";
    }
    return oss.str();
}
