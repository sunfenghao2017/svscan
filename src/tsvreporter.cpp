#include "stats.h"
#include "svrec.h"
#include "fusionopt.h"

void Stats::reportSVTSV(SVSet& svs, GeneInfoList& gl){
    std::ofstream fw(mOpt->tsvOut);
    fw << SVRec::gethead(mOpt->rnamode);
    for(uint32_t i = 0; i < gl.size(); ++i){
        // skip false positive insertions
        if(svs[i]->mSVT == 4 && mJctCnts[i].mFPIns > mSpnCnts[i].getAltDep() * mOpt->filterOpt->mMaxFPIns) continue;
        SVRec svr;
        // svType
        svr.svType = svutil::addID(svs[i]->mSVT);
        // svSize
        if(svs[i]->mSVT >= 5) svr.svSize = -1;
        else if(svs[i]->mSVT == 4) svr.svSize = svs[i]->mInsSeq.size();
        else svr.svSize = svs[i]->mSize;
        // bpMark
        svr.bpMark = svutil::getBpMark(svs[i]->mSVT);
        // bp1Chr bp1Pos bp2Chr bp2Pos
        svr.bp1Chr = gl[i].mChr1;
        svr.bp1Pos = gl[i].mPos1;
        svr.bp2Chr = gl[i].mChr2;
        svr.bp2Pos = gl[i].mPos2;
        // srCount dpCount srRescued dpRescued molRescued
        svr.srCount = svs[i]->mSRSupport;
        svr.dpCount = svs[i]->mPESupport;
        svr.srRescued = mJctCnts[i].getAltDep();
        svr.dpRescued = mSpnCnts[i].getAltDep();
        svr.molRescued = mTotalAltCnts[i];
        // srsrescued srsmalncnt srsmrate
        svr.srsrescued = svs[i]->mSRSResAllCnt;
        svr.srsmalncnt = svs[i]->mSRSResMAlnCnt;
        if(svr.srRescued > 0) svr.srsmrate = (double)(svr.srsmalncnt)/double(svr.srsrescued);
        else svr.srsmrate = 0;
        // srRefCount dpRefCount
        svr.srRefCount = mJctCnts[i].getRefDep();
        svr.dpRefCount = mSpnCnts[i].getRefDep();
        // AF
        if(svr.molRescued + std::max(svr.srRefCount, svr.dpRefCount)){
            svr.af = (double)(svr.molRescued)/(double)(svr.molRescued + std::max(svr.srRefCount, svr.dpRefCount));
        }else{
            svr.af = 0;
        }
        // insBp insSeq
        svr.insBp = svs[i]->mBpInsSeq.length();
        svr.insSeq = (svs[i]->mBpInsSeq.length() == 0 ? "-" : svs[i]->mBpInsSeq);
        // svSeq seqBp
        if(svs[i]->mSVT == 4){
            svr.svSeq = svs[i]->mInsSeq;
            svr.seqBp = 0;
        }else if(svs[i]->mPrecise){
            svr.svSeq = svs[i]->mConsensus;
            svr.seqBp = svs[i]->mGapCoord[0];
        }else{
            svr.svSeq = "-";
            svr.seqBp = 0;
        }
        // svID svtInt
        svr.id = svs[i]->mID;
        svr.svInt = svs[i]->mSVT;
        // bp1Gene bp2Gene
        if(gl[i].mGene1.empty()) svr.bp1Gene = "-";
        else svr.bp1Gene = gl[i].getTrs1();
        if(gl[i].mGene2.empty()) svr.bp2Gene = "-";
        else svr.bp2Gene = gl[i].getTrs2();
        // fuseGene fsMask fsHits
        svr.fuseGene = gl[i].getFuseGene();
        svr.fsMask = gl[i].getFsMask();
        svr.fsHits = svs[i]->mRealnRet;
        if(mOpt->rnamode){
            // ts1Name ts1Pos ts2Name ts2Pos
            svr.trs1Name = svs[i]->mNameChr1;
            svr.trs1Pos = svs[i]->mSVStart;
            svr.trs2Name = svs[i]->mNameChr2;
            svr.trs2Pos = svs[i]->mSVEnd;
            svr.rnamode = true;
        }
        // output
        fw << svr;
    }
    fw.close();
}

void Stats::makeFuseRec(const SVSet& svs, GeneInfoList& gl){
    mOpt->fuseOpt->init();
#ifdef DEBUG
    if(mOpt->debug & DEBUG_FOUTF){
        if(!mOpt->fuseOpt->mFsRptList.empty()){
            for(auto iter = mOpt->fuseOpt->mFusionRptMap.begin(); iter != mOpt->fuseOpt->mFusionRptMap.end(); ++iter){
                std::cout << iter->first << std::endl;
                std::cout << iter->second << std::endl;
            }
        }
    }
#endif
    // annotate extra gene fusion events
    if(!mOpt->fuseOpt->mExtraAnnoList.empty()){
        for(uint32_t i = 0; i < gl.size(); ++i){
            TrsRecList exgs = mOpt->fuseOpt->mExtraAnnotator.anno(svs[i]->mNameChr1, svs[i]->mSVStart, svs[i]->mSVStart + 1);
            TrsRecList exge = mOpt->fuseOpt->mExtraAnnotator.anno(svs[i]->mNameChr2, svs[i]->mSVEnd, svs[i]->mSVEnd + 1);
            for(uint32_t k = 0; k < exgs.size(); ++k){
                exgs[k].pos = svs[i]->mSVStart;
            }
            for(uint32_t k = 0; k < exge.size(); ++k){
                exge[k].pos = svs[i]->mSVEnd;
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
                    if(mOpt->fuseOpt->mExtraAnnotator.matchp(fg.hgene, fg.tgene)) extfl.push_back(fg);
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
                    if(mOpt->fuseOpt->mExtraAnnotator.matchp(fg.hgene, fg.tgene)) extfl.push_back(fg);
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
                gl[i].mFuseGene[j].status |= FUSION_FWITHMIRROR;
            }
        }
    }
    // mask other status of fusions
    for(uint32_t i = 0; i < gl.size(); ++i){
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
#ifdef DEBUG
            if(mOpt->debug & DEBUG_FOUTF){
                std::cout << gl[i].mFuseGene[j].debugStr();
            }
#endif
            if(mOpt->fuseOpt->hasBlackGene(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                gl[i].mFuseGene[j].status |= FUSION_FBLACKGENE;
            }
            if(mOpt->fuseOpt->inBlackList(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                gl[i].mFuseGene[j].status |= FUSION_FBLACKPAIR;
            }
            if(!mOpt->fuseOpt->validSV(svs[i]->mSVT, svs[i]->mNameChr1, svs[i]->mNameChr2, svs[i]->mSVStart, svs[i]->mSVEnd)){
                gl[i].mFuseGene[j].status |= FUSION_FFBG;
            }
            if(mOpt->fuseOpt->inWhiteList(gl[i].mFuseGene[j].hgene, gl[i].mFuseGene[j].tgene)){
                gl[i].mFuseGene[j].status |= FUSION_FINDB;
            }
            if(mOpt->fuseOpt->inWhiteList(gl[i].mFuseGene[j].tgene, gl[i].mFuseGene[j].hgene)){
                gl[i].mFuseGene[j].status |= FUSION_FMINDB;
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
                if(mOpt->fuseOpt->inSameSVRngMap(gname, exlist, svs[i]->mSVT)){
                    gl[i].mFuseGene[j].status &= (~(FUSION_FTOOSMALLSIZE | FUSION_FINSAMEGENE));
                }
            }
            gl[i].mFuseGene[j].status &= (~(FUSION_FERRREALN | FUSION_FMULTREALN));
            if(svs[i]->mPrecise){
                gl[i].mFuseGene[j].status |= FUSION_FPRECISE;
                if(svutil::simpleSeq(svs[i]->mConsensus.substr(0, svs[i]->mGapCoord[0])) ||
                   svutil::simpleSeq(svs[i]->mConsensus.substr(svs[i]->mGapCoord[1])) ||
                   svutil::tandemRepSeq(svs[i]->mConsensus, TandemRepeatThresholdMap)){
                    gl[i].mFuseGene[j].status |= FUSION_FLOWCOMPLEX;
                }
                if(svs[i]->mRealnRet < 0) gl[i].mFuseGene[j].status |= FUSION_FERRREALN;
                if(gl[i].mFuseGene[j].status & (FUSION_FINDB | FUSION_FMINDB)){
                    if(svs[i]->mRealnRet > mOpt->fuseOpt->mWhiteFilter.mMaxRepHit){
                        gl[i].mFuseGene[j].status |= FUSION_FMULTREALN;
                    }
                }else{
                    if(svs[i]->mRealnRet > mOpt->fuseOpt->mUsualFilter.mMaxRepHit){
                        gl[i].mFuseGene[j].status |= FUSION_FMULTREALN;
                    }
                }
            }
            if(!(gl[i].mFuseGene[j].status & (FUSION_FERRREALN | FUSION_FMULTREALN))){
                gl[i].mFuseGene[j].status |= FUSION_FPASSREALN;
            }
            std::string gene1, gene2;
            if(gl[i].mFuseGene[j].hfrom1){
                gene1 = gl[i].mGene1[gl[i].mFuseGene[j].hidx].gene;
                gene2 = gl[i].mGene2[gl[i].mFuseGene[j].tidx].gene;
            }else{
                gene1 = gl[i].mGene2[gl[i].mFuseGene[j].hidx].gene;
                gene2 = gl[i].mGene1[gl[i].mFuseGene[j].tidx].gene;
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
                bool omidb = (gl[i].mFuseGene[j].status & FUSION_FMINDB);
                bool oindb = (gl[i].mFuseGene[j].status & FUSION_FINDB);
                if(omidb ^ oindb){
                    if(omidb){
                        gl[i].mFuseGene[j].status &= (~FUSION_FMINDB);
                        gl[i].mFuseGene[j].status |= FUSION_FINDB;
                    }
                    if(oindb){
                        gl[i].mFuseGene[j].status &= (~FUSION_FINDB);
                        gl[i].mFuseGene[j].status |= FUSION_FMINDB;
                    }
                }
            }
        }
    }
    // fsCigar (RNA Only)
    if(mOpt->rnamode){
        for(uint32_t i = 0; i < gl.size(); ++i){
            for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
                if(gl[i].mFuseGene[j].tfrom1){
                    gl[i].mFuseGene[j].cigar = svutil::bp2cigar(gl[i].mGene2[gl[i].mFuseGene[j].hidx], gl[i].mGene1[gl[i].mFuseGene[j].tidx]);
                }else{
                    gl[i].mFuseGene[j].cigar = svutil::bp2cigar(gl[i].mGene1[gl[i].mFuseGene[j].hidx], gl[i].mGene2[gl[i].mFuseGene[j].tidx]);
                }
            }
            gl[i].mFsCigar = gl[i].mFuseGene[0].cigar;
        }
    }
    // construct fusionrecord
    for(uint32_t i = 0; i < gl.size(); ++i){
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            FusionRecord fsr;
            toFuseRec(fsr, svs[i], gl[i], j);
            fsr.maskFusion(mOpt->fuseOpt);
            gl[i].mFuseGene[j].status |= fsr.fsmask;
        }
    }
    // mask report mask
    for(uint32_t i = 0; i < gl.size(); ++i){
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            if((gl[i].mFuseGene[j].status & FUSION_FINDB) && (gl[i].mFuseGene[j].status & mOpt->fuseOpt->mIDBDropMask)) continue;
            if((!(gl[i].mFuseGene[j].status & FUSION_FINDB)) && (gl[i].mFuseGene[j].status & mOpt->fuseOpt->mNDBDropMask)) continue;
            for(auto& e: mOpt->fuseOpt->mKeepMasks){
                if((gl[i].mFuseGene[j].status & e) == e){
                    if((gl[i].mFuseGene[j].status & mOpt->fuseOpt->mPrimaryMask) == mOpt->fuseOpt->mPrimaryMask){
                        gl[i].mFuseGene[j].status |= FUSION_FPRIMARY;
                    }else{
                        gl[i].mFuseGene[j].status |= FUSION_FSUPPLEMENTARY;
                    }
                    break;
                }
            }
        }
    }
}

void Stats::reportFusionTSV(const SVSet& svs, GeneInfoList& gl){
    // get valid fusion list
    FusionRecordList frl;
    for(uint32_t i = 0; i < gl.size(); ++i){
        bool reported = false;
        for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
            if(gl[i].mFuseGene[j].status & FUSION_FPRIMARY){
                FusionRecord fsr;
                toFuseRec(fsr, svs[i], gl[i], j);
                frl.push_back(fsr);
                reported = true;
            }
        }
        if(!reported){
            for(uint32_t j = 0; j < gl[i].mFuseGene.size(); ++j){
                if(gl[i].mFuseGene[j].status & FUSION_FSUPPLEMENTARY){
                    FusionRecord fsr;
                    toFuseRec(fsr, svs[i], gl[i], j);
                    frl.push_back(fsr);
                }
            }
        }
    }
    std::sort(frl.begin(), frl.end());
    // mark report 
    uint32_t l = 0;
    for(uint32_t i = 1; i < frl.size(); ++i){
        if(!frl[i].samefs(frl[l])) frl[l].report = true;
        l = i;
    }
    if(!frl.empty()) frl[frl.size() - 1].report = true;
#ifdef DEBUG
    if(mOpt->debug & DEBUG_FOUTF){
        std::cout << "DEBUG_SORTED_FUSION_EVENTS:" << std::endl;
        for(uint32_t i = 0; i < frl.size(); ++i){
            std::cout << std::boolalpha << frl[i].report << "\t" <<  frl[i] << std::endl;
        }
    }
#endif
    // mark mirror fusion to keep optimal one in output
    markMirrorFusionEvent(frl);
    // output valid fusions
    std::string header = FusionRecord::gethead(mOpt->rnamode);
    std::ofstream fw(mOpt->fuseOpt->mOutFile);
    fw << header;
    for(uint32_t i = 0; i < frl.size(); ++i){
        if(frl[i].report) fw << frl[i];
    }
    fw.close();
}

void Stats::toFuseRec(FusionRecord& fsr, const SVRecord* svr, GeneInfo& gi, int32_t i){
    std::stringstream oss;
    fsr.fusionreads = mTotalAltCnts[svr->mID];
    fsr.totalreads = std::max(mJctCnts[svr->mID].getRefDep(), mSpnCnts[svr->mID].getRefDep()) + fsr.fusionreads;
    fsr.fuserate = (double)(fsr.fusionreads)/(double)(fsr.totalreads);
    fsr.fusegene = gi.mFuseGene[i].hgene + "->" + gi.mFuseGene[i].tgene; // FusionGene
    // FusionPattern
    if(gi.mFuseGene[i].hfrom1) fsr.fusepattern += gi.mGene1[gi.mFuseGene[i].hidx].strand;
    else fsr.fusepattern += gi.mGene2[gi.mFuseGene[i].hidx].strand;
    if(gi.mFuseGene[i].tfrom1) fsr.fusepattern += gi.mGene1[gi.mFuseGene[i].tidx].strand;
    else fsr.fusepattern += gi.mGene2[gi.mFuseGene[i].tidx].strand;
    // Gene1 Chr1 JunctionPosition1 Strand1 Transcript1
    if(gi.mFuseGene[i].hfrom1){
        fsr.gene1 = gi.mGene1[gi.mFuseGene[i].hidx].gene;
        fsr.chr1 = gi.mGene1[gi.mFuseGene[i].hidx].chr;
        fsr.jctpos1 = gi.mGene1[gi.mFuseGene[i].hidx].pos;
        fsr.strand1 = gi.mGene1[gi.mFuseGene[i].hidx].strand;
        fsr.transcript1 = gi.mGene1[gi.mFuseGene[i].hidx].getTrs();
        fsr.exon1 = gi.mGene1[gi.mFuseGene[i].hidx].exon;
    }else{
        fsr.gene1 = gi.mGene2[gi.mFuseGene[i].hidx].gene;
        fsr.chr1 = gi.mGene2[gi.mFuseGene[i].hidx].chr;
        fsr.jctpos1 = gi.mGene2[gi.mFuseGene[i].hidx].pos;
        fsr.strand1 = gi.mGene2[gi.mFuseGene[i].hidx].strand;
        fsr.transcript1 = gi.mGene2[gi.mFuseGene[i].hidx].getTrs();
        fsr.exon1 = gi.mGene2[gi.mFuseGene[i].hidx].exon;
    }
    // Gene2 Chr2 JunctionPosition2 Strand2 Transcript2
    if(gi.mFuseGene[i].tfrom1){
        fsr.gene2 = gi.mGene1[gi.mFuseGene[i].tidx].gene;
        fsr.chr2 = gi.mGene1[gi.mFuseGene[i].tidx].chr;
        fsr.jctpos2 = gi.mGene1[gi.mFuseGene[i].tidx].pos;
        fsr.strand2 = gi.mGene1[gi.mFuseGene[i].tidx].strand;
        fsr.transcript2 = gi.mGene1[gi.mFuseGene[i].tidx].getTrs();
        fsr.exon2 = gi.mGene1[gi.mFuseGene[i].tidx].exon;
    }else{
        fsr.gene2 = gi.mGene2[gi.mFuseGene[i].tidx].gene;
        fsr.chr2 = gi.mGene2[gi.mFuseGene[i].tidx].chr;
        fsr.jctpos2 = gi.mGene2[gi.mFuseGene[i].tidx].pos;
        fsr.strand2 = gi.mGene2[gi.mFuseGene[i].tidx].strand;
        fsr.transcript2 = gi.mGene2[gi.mFuseGene[i].tidx].getTrs();
        fsr.exon2 = gi.mGene2[gi.mFuseGene[i].tidx].exon;
    }
    // FusinSequence fseqBp
    if(svr->mSVT == 4 || (!svr->mPrecise)){
        fsr.fusionsequence = "-";
        fsr.fseqbp = 0;
    }else{
        fsr.fusionsequence = svr->mConsensus;
        fsr.fseqbp = svr->mGapCoord[0];
    }
    if(gi.mFuseGene[i].status & (FUSION_FINDB | FUSION_FMINDB)) fsr.indb = "Y"; // inDB
    else fsr.indb = "N";
    fsr.svt = svutil::addID(svr->mSVT);                       // svType
    if(svr->mSVT >= 5) fsr.svsize = -1;                       // svSize
    else fsr.svsize = svr->mSize;
    fsr.srcount = svr->mSRSupport;                            // srCount
    fsr.dpcount = svr->mPESupport;                            // dpCount
    fsr.srrescued =  mJctCnts[svr->mID].getAltDep();          // srRescued
    fsr.dprescued = mSpnCnts[svr->mID].getAltDep();           // dpRescued
    fsr.srrefcount = mJctCnts[svr->mID].getRefDep();          // srRefCount
    fsr.dprefcount = mSpnCnts[svr->mID].getRefDep();          // dpRefCount
    fsr.srsrescued = svr->mSRSResAllCnt;                      // srSRescued
    fsr.srsmalncnt = svr->mSRSResMAlnCnt;                     // srSResMaln
    if(fsr.srsrescued > 0){                                    // srSResMalnRate
        fsr.srsmrate = (double)(fsr.srsmalncnt)/fsr.srsrescued;
    }else{
        fsr.srsmrate = 0;
    }
    fsr.insbp = svr->mBpInsSeq.length();                      // insBp
    if(svr->mBpInsSeq.empty()) fsr.insseq = "-";              // insSeq
    else fsr.insseq = svr->mBpInsSeq;
    fsr.svid = svr->mID;                                      // svID
    fsr.svint = svr->mSVT;                                    // svtInt
    fsr.fsmask = gi.mFuseGene[i].status;                      // fsMask
    fsr.fsHits = svr->mRealnRet;                              // fsHits
    if(mOpt->rnamode){
        fsr.ts1name = svr->mNameChr1;                         // ts1Name
        fsr.ts1pos = svr->mSVStart;                           // ts1Pos
        fsr.ts2name = svr->mNameChr2;                         // ts2Name
        fsr.ts2pos = svr->mSVEnd;                             // ts2Pos
        fsr.cigar = gi.mFsCigar;                              // fsCigar
        if(gi.mGene1[i].gene != gi.mFuseGene[i].hgene){
            std::swap(fsr.ts1name, fsr.ts2name);
            std::swap(fsr.ts1pos, fsr.ts2pos);
        }
    }
    fsr.distance = mOpt->fuseOpt->geneNear(fsr.gene1, fsr.chr1, fsr.jctpos1, fsr.gene2, fsr.chr2); // distance
}
