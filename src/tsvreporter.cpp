#include "stats.h"
#include "svrec.h"
#include "fusionopt.h"

void Stats::reportSVTSV(SVSet& svs, GeneInfoList& gl){
    std::ofstream fw(mOpt->tsvOut);
    fw << SVRec::gethead(mOpt->rnamode);
    for(uint32_t i = 0; i < gl.size(); ++i){
        // skip false positive insertions
        if(svs[i].mSVT == 4 && mJctCnts[i].mFPIns > mSpnCnts[i].getAltDep() * mOpt->filterOpt->mMaxFPIns) continue;
        SVRec svr;
        // svType
        svr.svType = svutil::addID(svs[i].mSVT);
        // svSize
        if(svs[i].mSVT >= 5) svr.svSize = -1;
        else if(svs[i].mSVT == 4) svr.svSize = svs[i].mInsSeq.size();
        else svr.svSize = svs[i].mSize;
        // bpMark
        svr.bpMark = svutil::getBpMark(svs[i].mSVT);
        // bp1Chr bp1Pos bp2Chr bp2Pos
        svr.bp1Chr = gl[i].mChr1;
        svr.bp1Pos = gl[i].mPos1;
        svr.bp2Chr = gl[i].mChr2;
        svr.bp2Pos = gl[i].mPos2;
        // srCount dpCount srRescued dpRescued
        svr.srCount = svs[i].mSRSupport;
        svr.dpCount = svs[i].mPESupport;
        svr.srRescued = mJctCnts[i].getAltDep();
        svr.dpRescued = mSpnCnts[i].getAltDep();
        // srRefCount dpRefCount
        svr.srRefCount = mJctCnts[i].getRefDep();
        svr.dpRefCount = mSpnCnts[i].getRefDep();
        // AF
        if(svs[i].mPrecise){
            int32_t srv = std::max(mJctCnts[i].getAltDep(), svs[i].mSRSupport);
            svr.af = (double)(srv)/(mJctCnts[i].getRefDep() + srv);
        }else{
            int32_t dpv = std::max(mSpnCnts[i].getAltDep(), svs[i].mPESupport);
            svr.af = (double)(dpv)/(mSpnCnts[i].getRefDep() + dpv);
        }
        // insBp insSeq
        svr.insBp = svs[i].mBpInsSeq.length();
        svr.insSeq = (svs[i].mBpInsSeq.length() == 0 ? "-" : svs[i].mBpInsSeq);
        // svSeq seqBp
        if(svs[i].mSVT == 4){
            svr.svSeq = svs[i].mInsSeq;
            svr.seqBp = 0;
        }else if(svs[i].mPrecise){
            svr.svSeq = svs[i].mConsensus;
            svr.seqBp = svs[i].mGapCoord[0];
        }else{
            svr.svSeq = "-";
            svr.seqBp = 0;
        }
        // svID svtInt
        svr.id = svs[i].mID;
        svr.svInt = svs[i].mSVT;
        // bp1Gene bp2Gene
        if(gl[i].mGene1.empty()) svr.bp1Gene = "-";
        else svr.bp1Gene = gl[i].getTrs1();
        if(gl[i].mGene2.empty()) svr.bp2Gene = "-";
        else svr.bp2Gene = gl[i].getTrs2();
        // fuseGene fsMask
        svr.fuseGene = gl[i].getFuseGene();
        svr.fsMask = gl[i].getFsMask();
        if(mOpt->rnamode){
            // ts1Name ts1Pos ts2Name ts2Pos
            svr.trs1Name = svs[i].mNameChr1;
            svr.trs1Pos = svs[i].mSVStart;
            svr.trs2Name = svs[i].mNameChr2;
            svr.trs2Pos = svs[i].mSVEnd;
            svr.rnamode = true;
        }
        // output
        fw << svr;
    }
    fw.close();
}

void Stats::maskFuseRec(const SVSet& svs, GeneInfoList& gl){
    mOpt->fuseOpt->init();
    RealnFilter realnFilter(mOpt->alnref);
    // annotate extra gene fusion events
    if(!mOpt->fuseOpt->mExtraAnnoList.empty()){
        for(uint32_t i = 0; i < gl.size(); ++i){
            TrsRecList exgs = mOpt->fuseOpt->mExtraAnnotator.anno(svs[i].mNameChr1, svs[i].mSVStart, svs[i].mSVStart + 1);
            TrsRecList exge = mOpt->fuseOpt->mExtraAnnotator.anno(svs[i].mNameChr2, svs[i].mSVEnd, svs[i].mSVEnd + 1);
            for(uint32_t k = 0; k < exgs.size(); ++k){
                exgs[k].pos = svs[i].mSVStart;
            }
            for(uint32_t k = 0; k < exge.size(); ++k){
                exge[k].pos = svs[i].mSVEnd;
            }
            int32_t g1osize = gl[i].mGene1.size();
            int32_t g2osize = gl[i].mGene2.size();
            std::copy(exgs.begin(), exgs.end(), std::back_inserter(gl[i].mGene1));
            std::copy(exge.begin(), exge.end(), std::back_inserter(gl[i].mGene2));
            FuseGeneList extfl;
            for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
                for(uint32_t k = 0; k < exgs.size(); ++k){
                    FuseGene fg = gl[i].mFuseGene[j];
                    if(fg.hfrom1){
                        fg.hgene = exgs[k].gene;
                        fg.hidx = g1osize + k;
                    }else{
                        fg.tgene = exgs[k].gene;
                        fg.tidx = g1osize + k;
                    }
                    extfl.push_back(fg);
                }
                for(uint32_t k = 0; k < exge.size(); ++k){
                    FuseGene fg = gl[i].mFuseGene[j];
                    if(fg.tfrom1){
                        fg.hgene = exge[k].gene;
                        fg.hidx = g2osize + k;
                    }else{
                        fg.tgene = exge[k].gene;
                        fg.tidx = g2osize + k;
                    }
                    extfl.push_back(fg);
                }
            }
            for(uint32_t k = 0; k < extfl.size(); ++k){
                if(extfl[k].hgene != "-" && extfl[k].tgene != "-"){
                    extfl[k].status |= FUSION_FALLGENE;
                    extfl[k].status |= FUSION_FCOMMONHOTDIRECT;
                    extfl[k].status |= FUSION_FINDB;
                    extfl[k].status |= FUSION_FNORMALCATDIRECT;
                }
            }
            std::copy(extfl.begin(), extfl.end(), std::back_inserter(gl[i].mFuseGene));
        }
    }
    // mask hot gene status
    std::map<std::string, std::set<std::string>> fpairs;
    for(uint32_t i = 0; i < gl.size(); ++i){
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            if(mOpt->rnamode) gl[i].mFuseGene[j].status |= FUSION_FCALLFROMRNASEQ; // mask rna/dna calling
            if(svs[i].mPassRealn) gl[i].mFuseGene[j].status |= FUSION_FREALNPASSED;
            if(mOpt->fuseOpt->hasWhiteGene(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                gl[i].mFuseGene[j].status |= FUSION_FHOTGENE;
                if(mOpt->fuseOpt->matchHotDirec(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                    gl[i].mFuseGene[j].status |= FUSION_FCOMMONHOTDIRECT;
                }
            }
            if(gl[i].mFuseGene[j].hgene != gl[i].mFuseGene[j].tgene){
                fpairs[gl[i].mFuseGene[j].hgene].insert(gl[i].mFuseGene[j].tgene);
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
            if(mOpt->debug & DEBUG_FOUTF){
                std::cout << gl[i].mFuseGene[j].debugStr();
            }
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
                int32_t minExon = -1, maxExon = -2;
                std::vector<int32_t> exlist;
                std::string gname;
                if(gl[i].mFuseGene[j].hfrom1){
                    if(gl[i].mGene1[gl[i].mFuseGene[j].hidx].near(gl[i].mGene2[gl[i].mFuseGene[j].tidx], mOpt->filterOpt->mMinDelRatio)){
                        gl[i].mFuseGene[j].status |= FUSION_FTOOSMALLSIZE;
                    }
                    minExon = std::min(std::atoi(gl[i].mGene1[gl[i].mFuseGene[j].hidx].number.c_str()), 
                                       std::atoi(gl[i].mGene2[gl[i].mFuseGene[j].tidx].number.c_str()));
                    maxExon = std::max(std::atoi(gl[i].mGene1[gl[i].mFuseGene[j].hidx].number.c_str()),
                                       std::atoi(gl[i].mGene2[gl[i].mFuseGene[j].tidx].number.c_str()));
                    gname = gl[i].mGene1[gl[i].mFuseGene[j].hidx].gene;
                }else{
                    if(gl[i].mGene1[gl[i].mFuseGene[j].tidx].near(gl[i].mGene2[gl[i].mFuseGene[j].hidx], mOpt->filterOpt->mMinDelRatio)){
                        gl[i].mFuseGene[j].status |= FUSION_FTOOSMALLSIZE;
                    }
                    minExon = std::min(std::atoi(gl[i].mGene1[gl[i].mFuseGene[j].tidx].number.c_str()), 
                                       std::atoi(gl[i].mGene2[gl[i].mFuseGene[j].hidx].number.c_str()));
                    maxExon = std::max(std::atoi(gl[i].mGene1[gl[i].mFuseGene[j].tidx].number.c_str()),
                                       std::atoi(gl[i].mGene2[gl[i].mFuseGene[j].hidx].number.c_str()));
                    gname = gl[i].mGene1[gl[i].mFuseGene[j].tidx].gene;
                }
                for(int32_t ecnt = minExon; ecnt <= maxExon; ++ecnt) exlist.push_back(ecnt);
                if(!mOpt->fuseOpt->inSameSVRngMap(gname, exlist, svs[i].mSVT)){
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
            int32_t srv = std::max((int32_t)mJctCnts[i].getAltDep(), svs[i].mSRSupport);
            int32_t srr = mJctCnts[i].getRefDep();
            int32_t dpv = std::max((int32_t)mSpnCnts[i].getAltDep(), svs[i].mPESupport);
            int32_t dpr = mSpnCnts[i].getRefDep();
            if(svs[i].mPrecise) af = (double)(srv)/(double)(srv + srr);
            else af = (double)(dpv)/(double)(dpv + dpr);
            if(gl[i].mFuseGene[j].status & (FUSION_FINDB | FUSION_FMIRRORINDB)){// fusion in public database
                if(svs[i].mPrecise){
                    if(srv < mOpt->fuseOpt->mWhiteFilter.mMinSupport && dpv < mOpt->fuseOpt->mWhiteFilter.mMinSupport){
                        gl[i].mFuseGene[j].status |= FUSION_FLOWSUPPORT;
                    }
                    if(svs[i].mSRSupport < mOpt->fuseOpt->mWhiteFilter.mMinSRSeed && svs[i].mPESupport < mOpt->fuseOpt->mWhiteFilter.mMinDPSeed){
                        gl[i].mFuseGene[j].status |= FUSION_FLOWSUPPORT;
                    }
                }else{
                    if(dpv < mOpt->fuseOpt->mWhiteFilter.mMinSupport){
                        gl[i].mFuseGene[j].status |= FUSION_FLOWSUPPORT;
                    }
                    if(svs[i].mPESupport < mOpt->fuseOpt->mWhiteFilter.mMinDPSeed){
                        gl[i].mFuseGene[j].status |= FUSION_FLOWSUPPORT;
                    }
                }
                if(af < mOpt->fuseOpt->mWhiteFilter.mMinVAF){
                    gl[i].mFuseGene[j].status |= FUSION_FLOWAF;
                }
                if(svs[i].mPrecise){
                    if((srv + srr) < mOpt->fuseOpt->mWhiteFilter.mMinDepth && (dpv + dpr) < mOpt->fuseOpt->mWhiteFilter.mMinDepth){
                        gl[i].mFuseGene[j].status |= FUSION_FLOWDEPTH;
                    }
                }else{
                    if((dpr + dpv) < mOpt->fuseOpt->mWhiteFilter.mMinDepth){
                        gl[i].mFuseGene[j].status |= FUSION_FLOWDEPTH;
                    }
                }
                if((svs[i].mSVT != 4) && gl[i].mFuseGene[j].status & FUSION_FINSAMEGENE){
                    if(svs[i].mSize < mOpt->fuseOpt->mWhiteFilter.mMinIntraGeneSVSize){
                        gl[i].mFuseGene[j].status |= FUSION_FTOOSMALLSIZE;
                    }
                }
            }else if(gl[i].mFuseGene[j].status & FUSION_FHOTGENE){// fusion in whitelist
                if(svs[i].mPrecise){
                    if(srv < mOpt->fuseOpt->mUsualFilter.mMinSupport && dpv < mOpt->fuseOpt->mUsualFilter.mMinSupport){
                        gl[i].mFuseGene[j].status |= FUSION_FLOWSUPPORT;
                    }
                    if(svs[i].mSRSupport < mOpt->fuseOpt->mUsualFilter.mMinSRSeed && svs[i].mPESupport < mOpt->fuseOpt->mUsualFilter.mMinDPSeed){
                        gl[i].mFuseGene[j].status |= FUSION_FLOWSUPPORT;
                    }
                }else{
                    gl[i].mFuseGene[j].status |= FUSION_FLOWSUPPORT;
                }
                if(af < mOpt->fuseOpt->mUsualFilter.mMinVAF){
                    gl[i].mFuseGene[j].status |= FUSION_FLOWAF;
                }
                if(svs[i].mPrecise){
                    if((srv + srr) < mOpt->fuseOpt->mUsualFilter.mMinDepth && (dpv + dpr) < mOpt->fuseOpt->mUsualFilter.mMinDepth){
                        gl[i].mFuseGene[j].status |= FUSION_FLOWDEPTH;
                    }
                }else{
                    gl[i].mFuseGene[j].status |= FUSION_FLOWDEPTH;
                }
                if((svs[i].mSVT != 4) && gl[i].mFuseGene[j].status & FUSION_FINSAMEGENE){
                    if(svs[i].mSize < mOpt->fuseOpt->mUsualFilter.mMinIntraGeneSVSize){
                        gl[i].mFuseGene[j].status |= FUSION_FTOOSMALLSIZE;
                    }
                }
            }
        }   
    }
    // adjust hot gene in normal catenation direction
    for(uint32_t i = 0; i < gl.size(); ++i){
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            if(gl[i].mFuseGene[j].status & FUSION_FNORMALCATDIRECT) continue;
            if(gl[i].mFuseGene[j].status & FUSION_FCOMMONHOTDIRECT) continue;
            if(gl[i].mFuseGene[j].status & FUSION_FHOTGENE){
                std::swap(gl[i].mFuseGene[j].hend, gl[i].mFuseGene[j].tend);
                std::swap(gl[i].mFuseGene[j].hfrom1, gl[i].mFuseGene[j].tfrom1);
                std::swap(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene);
                std::swap(gl[i].mFuseGene[j].hidx, gl[i].mFuseGene[j].tidx);
                std::swap(gl[i].mFuseGene[j].hstrand, gl[i].mFuseGene[j].tstrand);
                gl[i].mFuseGene[j].status |= FUSION_FCOMMONHOTDIRECT;
                if(gl[i].mFuseGene[j].status & FUSION_FHTFLSWAPPED){
                    gl[i].mFuseGene[j].status &= (~FUSION_FHTFLSWAPPED);
                }else{
                    gl[i].mFuseGene[j].status |= FUSION_FHTFLSWAPPED;
                }
            }
        }
    }
    // drop bits mask of all fusion events, if an fusion match any bit in FUSION_DROP_MASK, it will not be reported
    TFUSION_FLAG FUSION_DROP_MASK = (FUSION_FBLACKGENE | FUSION_FBLACKPAIR  | FUSION_FFBG | FUSION_FLOWCOMPLEX |
                                     FUSION_FTOOSMALLSIZE | FUSION_FLOWAF | FUSION_FLOWSUPPORT | FUSION_FLOWDEPTH);
    // primary keep bits mask, fusion reported as primary must match all the bits in PRIMARY_KEEP_MASK
    TFUSION_FLAG PRIMARY_KEEP_MASK = (FUSION_FNORMALCATDIRECT | FUSION_FCOMMONHOTDIRECT | FUSION_FINDB);
    // keep bits mask, an fusion to be reported must match all bits in FUSION_KEEP_MASK
    TFUSION_FLAG FUSION_KEEP_MASK = (FUSION_FHOTGENE | FUSION_FREALNPASSED | FUSION_FALLGENE);
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
    std::string header = FusionRecord::gethead(mOpt->rnamode);
    std::ofstream fw(mOpt->fuseOpt->mOutFile);
    std::ofstream fs(mOpt->fuseOpt->mSupFile);
    fw << header;
    fs << header;
    for(uint32_t i = 0; i < gl.size(); ++i){
        bool reported = false;
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            if(gl[i].mFuseGene[j].status & FUSION_FPRIMARY){
                FusionRecord fsr;
                toFuseRec(fsr, svs[i], gl[i], j);
                fw << fsr;
                reported = true;
            }
        }
        if(!reported){
            for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
                if(gl[i].mFuseGene[j].status & FUSION_FSUPPLEMENTARY){
                    FusionRecord fsr;
                    toFuseRec(fsr, svs[i], gl[i], j);
                    fs << fsr;
                }
            }
        }
    }
    fw.close();
    fs.close();
}

void Stats::toFuseRec(FusionRecord& fsr, SVRecord& svr, GeneInfo& gi, int32_t i){
    std::stringstream oss;
    float af = 0.0;
    int32_t srv = std::max((int32_t)mJctCnts[svr.mID].getAltDep(), svr.mSRSupport);
    int32_t srr = mJctCnts[svr.mID].getRefDep();
    int32_t dpv = mSpnCnts[svr.mID].getAltDep();
    int32_t dpr = mSpnCnts[svr.mID].getRefDep();
    if(svr.mPrecise) af = (double)(srv)/(double)(srv + srr);
    else af = (double)(dpv)/(double)(dpv + dpr);
    fsr.fusegene = gi.mFuseGene[i].hgene + "->" + gi.mFuseGene[i].tgene; // FusionGene
    // FusionPattern
    if(gi.mFuseGene[i].hfrom1) fsr.fusepattern += gi.mGene1[gi.mFuseGene[i].hidx].strand;
    else fsr.fusepattern += gi.mGene2[gi.mFuseGene[i].hidx].strand;
    if(gi.mFuseGene[i].tfrom1) fsr.fusepattern += gi.mGene1[gi.mFuseGene[i].tidx].strand;
    else fsr.fusepattern += gi.mGene2[gi.mFuseGene[i].tidx].strand;
    if(svr.mPrecise){// FusionReads TotalReads FusionRate
        fsr.fusionreads = srv;
        fsr.totalreads = srv + srr;
        fsr.fuserate = af;
    }else{
        fsr.fusionreads = dpv;
        fsr.totalreads = dpr + dpv;
        fsr.fuserate = af;
    }
    // Gene1 Chr1 JunctionPosition1 Strand1 Transcript1
    if(gi.mFuseGene[i].hfrom1){
        fsr.gene1 = gi.mGene1[gi.mFuseGene[i].hidx].gene;
        fsr.chr1 = gi.mGene1[gi.mFuseGene[i].hidx].chr;
        fsr.junctionposition1 = gi.mGene1[gi.mFuseGene[i].hidx].pos;
        fsr.strand1 = gi.mGene1[gi.mFuseGene[i].hidx].strand;
        fsr.transcript1 = gi.mGene1[gi.mFuseGene[i].hidx].getTrs();
        fsr.exon1 = gi.mGene1[gi.mFuseGene[i].hidx].exon;
    }else{
        fsr.gene1 = gi.mGene2[gi.mFuseGene[i].hidx].gene;
        fsr.chr1 = gi.mGene2[gi.mFuseGene[i].hidx].chr;
        fsr.junctionposition1 = gi.mGene2[gi.mFuseGene[i].hidx].pos;
        fsr.strand1 = gi.mGene2[gi.mFuseGene[i].hidx].strand;
        fsr.transcript1 = gi.mGene2[gi.mFuseGene[i].hidx].getTrs();
        fsr.exon1 = gi.mGene2[gi.mFuseGene[i].hidx].exon;
    }
    // Gene2 Chr2 JunctionPosition2 Strand2 Transcript2
    if(gi.mFuseGene[i].tfrom1){
        fsr.gene2 = gi.mGene1[gi.mFuseGene[i].tidx].gene;
        fsr.chr2 = gi.mGene1[gi.mFuseGene[i].tidx].chr;
        fsr.junctionposition2 = gi.mGene1[gi.mFuseGene[i].tidx].pos;
        fsr.strand2 = gi.mGene1[gi.mFuseGene[i].tidx].strand;
        fsr.transcript2 = gi.mGene1[gi.mFuseGene[i].tidx].getTrs();
        fsr.exon2 = gi.mGene1[gi.mFuseGene[i].tidx].exon;
    }else{
        fsr.gene2 = gi.mGene2[gi.mFuseGene[i].tidx].gene;
        fsr.chr2 = gi.mGene2[gi.mFuseGene[i].tidx].chr;
        fsr.junctionposition2 = gi.mGene2[gi.mFuseGene[i].tidx].pos;
        fsr.strand2 = gi.mGene2[gi.mFuseGene[i].tidx].strand;
        fsr.transcript2 = gi.mGene2[gi.mFuseGene[i].tidx].getTrs();
        fsr.exon2 = gi.mGene2[gi.mFuseGene[i].tidx].exon;
    }
    // FusinSequence fseqBp
    if(svr.mSVT == 4 || (!svr.mPrecise)){
        fsr.fusionsequence = "-";
        fsr.fseqbp = 0;
    }else{
        fsr.fusionsequence = svr.mConsensus;
        fsr.fseqbp = svr.mGapCoord[0];
    }
    if(gi.mFuseGene[i].status & FUSION_FINDB) fsr.indb = "Y"; // inDB
    else fsr.indb = "N";
    fsr.svt = svutil::addID(svr.mSVT);                        // svType
    if(svr.mSVT >= 5) fsr.svsize = -1;                        // svSize
    else fsr.svsize = svr.mSize;
    fsr.srcount = svr.mSRSupport;                             // srCount
    fsr.dpcount = svr.mPESupport;                             // dpCount
    fsr.srrescued =  srv;                                     // srRescued
    fsr.dprescued = mSpnCnts[svr.mID].getAltDep();            // dpRescued
    fsr.srrefcount = mJctCnts[svr.mID].getRefDep();           // srRefCount
    fsr.dprefcount = mSpnCnts[svr.mID].getRefDep();           // dpRefCount
    fsr.insbp = svr.mBpInsSeq.length();                       // insBp
    if(svr.mBpInsSeq.empty()) fsr.insseq = "-";               // insSeq
    else fsr.insseq = svr.mBpInsSeq;
    fsr.svid = svr.mID;                                       // svID
    fsr.svint = svr.mSVT;                                     // svtInt
    fsr.fsmask = gi.mFuseGene[i].status;                      // fsMask
    if(mOpt->rnamode){
        fsr.ts1name = svr.mNameChr1;                          // ts1Name
        fsr.ts1pos = svr.mSVStart;                            // ts1Pos
        fsr.ts2name = svr.mNameChr2;                          // ts2Name
        fsr.ts2pos = svr.mSVEnd;                              // ts2Pos
        if(gi.mGene1[i].gene != gi.mFuseGene[i].hgene){
            std::swap(fsr.ts1name, fsr.ts2name);
            std::swap(fsr.ts1pos, fsr.ts2pos);
        }
        fsr.cigar = gi.mFuseGene[i].cigar;                    // fsCigar
    }
}
