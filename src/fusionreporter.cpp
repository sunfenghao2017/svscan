#include "fusionreporter.h"

FusionReporter::FusionReporter(){
    fuseOpt = new FusionOptions();
    softEnv = new Software();
    softEnv->cmp += "version: " + softEnv->ver + "\n";
    softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
    rnamode = false;
}

FusionReporter::~FusionReporter(){
    if(softEnv) delete softEnv;
    if(fuseOpt) delete fuseOpt;
}

void FusionReporter::update(int argc, char** argv){
    // update software environment records
    softEnv->cwd = util::cwd();
    for(int i = 0; i < argc; ++i){
        softEnv->cmd.append(argv[i]);
        softEnv->cmd.append(" ");
    }
}

void FusionReporter::str2trsl(TrsRecList& trsl, const std::string& trStr){
    if(trStr == "-") return;
    std::vector<std::string> vstr;
    util::split(trStr, vstr, ";");
    for(uint32_t i = 0; i < vstr.size(); ++i){
        std::vector<std::string> tvstr;
        util::split(vstr[i], tvstr, ",");
        TrsRec tr;
        tr.gene = tvstr[0];
        tr.name = tvstr[1];
        tr.unit = tvstr[2];
        tr.strand = tvstr[3];
        tr.number = tvstr[4];
        tr.exon = std::atoi(tvstr[5].c_str());
        trsl.push_back(tr);
    }
}

void FusionReporter::str2fsgs(FuseGeneList& fsgl, const std::string& fsStr, const std::string& fmsks, const TrsRecList& bp1trs, const TrsRecList& bp2trs){
    if(fsStr.empty()) return;
    std::vector<std::string> mvst;
    util::split(fmsks, mvst, ";");
    std::vector<std::string> vstr;
    util::split(fsStr, vstr, ";");
    for(uint32_t i = 0; i < vstr.size(); ++i){
        std::vector<std::string> fvstr;
        util::split(vstr[i], fvstr, "|");
        std::vector<std::string> hvstr, tvstr;
        util::split(fvstr[1], hvstr, ",");
        util::split(fvstr[2], tvstr, ",");
        FuseGene fg;
        fg.status = std::atoi(mvst[i].c_str());
        fg.hgene = hvstr[0];
        fg.hend = hvstr[1];
        fg.hstrand = hvstr[2];
        fg.tgene = tvstr[0];
        fg.tend = tvstr[1];
        fg.tstrand = tvstr[2];
        if(fg.status & FUSION_FHTFLSWAPPED){
            fg.hfrom1 = false;
            fg.tfrom1 = true;
        }else{
            fg.hfrom1 = true;
            fg.tfrom1 = false;
        }
        for(uint32_t b1 = 0; b1 < bp1trs.size(); ++b1){
            if(fg.hfrom1 && bp1trs[b1].gene == fg.hgene) fg.hidx = b1;
            if(fg.tfrom1 && bp1trs[b1].gene == fg.tgene) fg.tidx = b1;
        }
        for(uint32_t b2 = 0; b2 < bp2trs.size(); ++b2){
            if((!fg.hfrom1) && bp2trs[b2].gene == fg.hgene) fg.hidx = b2;
            if((!fg.tfrom1) && bp2trs[b2].gene == fg.tgene) fg.tidx = b2;
        }
        if(fvstr.size() > 2) fg.cigar = fvstr[2];
        fsgl.push_back(fg);
    }
}

void FusionReporter::report(){
    sv2fsl(fuseList);
    std::sort(fuseList.begin(), fuseList.end());
    // mark report
    uint32_t l = 0;
    for(uint32_t i = 1; i < fuseList.size(); ++i){
        if(!fuseList[i].samefs(fuseList[l])) fuseList[l].report = true;
        l = i;
    }
    if(!fuseList.empty()) fuseList[fuseList.size() - 1].report = true;
    // output valid fusions
    std::string header = FusionRecord::gethead(rnamode);
    std::ofstream fw(fuseOpt->mOutFile);
    std::ofstream fs(fuseOpt->mSupFile);
    fw << header;
    fs << header;
    for(auto& e: fuseList){
        if(e.report){
            if(e.fsmask & FUSION_FPRIMARY){
                fw << e;
            }
            if(e.fsmask & FUSION_FSUPPLEMENTARY){
                fs << e;
            }
        }
    }
    fw.close();
    fs.close();
}

void FusionReporter::sv2fsl(FusionRecordList& fsrl){
    fuseOpt->init();
    std::ofstream fsv;
    if(!fuseOpt->mSVModFile.empty()){ // update sv tsv file if needed
        fsv.open(fuseOpt->mSVModFile.c_str());
    }
    // drop bits mask of all fusion events, if an fusion match any bit in FUSION_DROP_MASK, it will not be reported
    TFUSION_FLAG FUSION_DROP_MASK = (FUSION_FBLACKGENE | FUSION_FBLACKPAIR  | FUSION_FFBG | FUSION_FLOWCOMPLEX |
                                     FUSION_FTOOSMALLSIZE | FUSION_FLOWAF | FUSION_FLOWSUPPORT | FUSION_FLOWDEPTH);
    // primary keep bits mask, fusion reported as primary must match all the bits in PRIMARY_KEEP_MASK
    TFUSION_FLAG PRIMARY_KEEP_MASK = (FUSION_FNORMALCATDIRECT | FUSION_FCOMMONHOTDIRECT | FUSION_FINDB | FUSION_FINREPORTRNG);
    // keep bits mask, an fusion to be reported must match all bits in FUSION_KEEP_MASK
    TFUSION_FLAG FUSION_KEEP_MASK = (FUSION_FHOTGENE | FUSION_FREALNPASSED | FUSION_FALLGENE);
    // supplementary fusion additional conditions
    std::ifstream fr(fuseOpt->mInfile);
    std::string tmpstr;
    std::vector<std::string> vstr;
    std::getline(fr, tmpstr);
    if(!fuseOpt->mSVModFile.empty()) fsv << tmpstr << "\n";
    while(std::getline(fr, tmpstr)){
        if(debug){
            std::cout << tmpstr << std::endl;
        }
        util::split(tmpstr, vstr, "\t");
        TrsRecList trsl1, trsl2;
        str2trsl(trsl1, vstr[20]);
        str2trsl(trsl2, vstr[21]);
        FuseGeneList fgl;
        str2fsgs(fgl, vstr[22], vstr[23], trsl1, trsl2);
        if(debug){
            std::cout << "FuseGeneList fgl: " << std::endl;
            for(uint32_t flidx = 0; flidx < fgl.size(); ++flidx){
                std::cout << fgl[flidx].debugStr() << std::endl;
            }
            std::cout << "TrsRecList trsl1: " << std::endl;
            for(uint32_t trsidx = 0; trsidx < trsl1.size(); ++trsidx){
                std::cout << trsl1[trsidx].toStr() << std::endl;
            }
            std::cout << "TrsRecList trsl2: " << std::endl;
            for(uint32_t trsidx = 0; trsidx < trsl2.size(); ++trsidx){
                std::cout << trsl2[trsidx].toStr() << std::endl;
            }
        }
        int32_t svt = std::atoi(vstr[19].c_str());
        std::string chr1 = vstr[3];
        std::string chr2 = vstr[5];
        int32_t start = std::atoi(vstr[4].c_str());
        int32_t end = std::atoi(vstr[6].c_str());
        int32_t srs = std::atoi(vstr[7].c_str());
        int32_t dps = std::atoi(vstr[8].c_str());
        int32_t srv = std::atoi(vstr[7].c_str());
        int32_t dpv = std::atoi(vstr[8].c_str());
        int32_t srr = std::atoi(vstr[11].c_str());
        int32_t dpr = std::atoi(vstr[12].c_str());
        float af = std::atof(vstr[13].c_str());
        bool notinbg= fuseOpt->validSV(svt, chr1, chr2, start, end);
        FusionRecordList frl;
        for(uint32_t i = 0; i < fgl.size(); ++i){
            FusionRecord fgr;
            fgr.fsmask = fgl[i].status;
            if(!rnamode && (fgr.fsmask & FUSION_FCALLFROMRNASEQ)) rnamode = true;
            fgr.cigar = fgl[i].cigar;
            fgr.gene1 = fgl[i].hgene;
            fgr.gene2 = fgl[i].tgene;
            if(svt == 4){
                fgl[i].hfrom1 = true;
                fgl[i].tfrom1 = false;
            }
            if(fgl[i].hfrom1) fgr.chr1 = chr1;
            else fgr.chr1 = chr2;
            if(fgl[i].tfrom1) fgr.chr2 = chr1;
            else fgr.chr2 = chr2;
            if(notinbg) fgr.fsmask &= (~FUSION_FFBG);
            else fgr.fsmask |= FUSION_FFBG;
            if(fuseOpt->hasWhiteGene(fgr.gene1, fgr.gene2)) fgr.fsmask |= FUSION_FHOTGENE;
            else fgr.fsmask &= (~FUSION_FHOTGENE);
            if(fuseOpt->inWhiteList(fgr.gene1, fgr.gene2)) fgr.fsmask |= FUSION_FINDB;
            else fgr.fsmask &= (~FUSION_FINDB);
            if(fuseOpt->inWhiteList(fgr.gene2, fgr.gene1)) fgr.fsmask |= FUSION_FMIRRORINDB;
            else fgr.fsmask &= (~FUSION_FMIRRORINDB);
            if(fuseOpt->hasBlackGene(fgr.gene1, fgr.gene2)) fgr.fsmask |= FUSION_FBLACKGENE;
            else fgr.fsmask &= (~FUSION_FBLACKGENE);
            if(fuseOpt->inBlackList(fgr.gene1, fgr.gene2)) fgr.fsmask |= FUSION_FBLACKPAIR;
            else fgr.fsmask &= (~FUSION_FBLACKPAIR);
            if(fuseOpt->matchHotDirec(fgr.gene1, fgr.gene2)) fgr.fsmask |= FUSION_FCOMMONHOTDIRECT;
            else fgr.fsmask &= (~FUSION_FCOMMONHOTDIRECT);
            if(fgl[i].status & FUSION_FINSAMEGENE && fgl[i].status & FUSION_FCALLFROMRNASEQ){
                std::vector<int32_t> exonl;
                int32_t minExon = -1, maxExon = -2;
                std::string gname;
                if(fgl[i].hfrom1){
                    minExon = std::min(std::atoi(trsl1[fgl[i].hidx].number.c_str()), std::atoi(trsl2[fgl[i].tidx].number.c_str()));
                    maxExon = std::max(std::atoi(trsl1[fgl[i].hidx].number.c_str()), std::atoi(trsl2[fgl[i].tidx].number.c_str()));
                    gname = trsl1[fgl[i].hidx].gene;
                }else{
                    minExon = std::min(std::atoi(trsl1[fgl[i].tidx].number.c_str()), std::atoi(trsl2[fgl[i].hidx].number.c_str()));
                    maxExon = std::max(std::atoi(trsl1[fgl[i].tidx].number.c_str()), std::atoi(trsl2[fgl[i].hidx].number.c_str()));
                    gname = trsl1[fgl[i].tidx].gene;
                }
                for(int32_t excnt = minExon; excnt <= maxExon; ++excnt) exonl.push_back(excnt);
                if(!fuseOpt->inSameSVRngMap(gname, exonl, fgr.svint)) fgr.fsmask |= FUSION_FTOOSMALLSIZE;
            }
            if(fgr.fsmask & (FUSION_FINDB | FUSION_FMIRRORINDB)){// fusion in public database
                if(srv){
                    if((srv < fuseOpt->mWhiteFilter.mMinSupport && dpv < fuseOpt->mWhiteFilter.mMinSupport) || 
                       (srs < fuseOpt->mWhiteFilter.mMinSRSeed && dps < fuseOpt->mWhiteFilter.mMinDPSeed)){
                        fgr.fsmask |= FUSION_FLOWSUPPORT;
                    }else{
                        fgr.fsmask &= (~FUSION_FLOWSUPPORT);
                    }
                }else{
                    if(dpv < fuseOpt->mWhiteFilter.mMinSupport || dps < fuseOpt->mWhiteFilter.mMinDPSeed){
                        fgr.fsmask |= FUSION_FLOWSUPPORT;
                    }else{
                        fgr.fsmask &= (~FUSION_FLOWSUPPORT);
                    }
                }
                if(af < fuseOpt->mWhiteFilter.mMinVAF){
                    fgr.fsmask |= FUSION_FLOWAF;
                }else{
                    fgr.fsmask &= (~FUSION_FLOWAF);
                }
                if(srv){
                    if((srv + srr) < fuseOpt->mWhiteFilter.mMinDepth){
                        fgr.fsmask |= FUSION_FLOWDEPTH;
                    }else{
                        fgr.fsmask &= (~FUSION_FLOWDEPTH);
                    }
                }else{
                    if((dpr + dpv) < fuseOpt->mWhiteFilter.mMinDepth){
                        fgr.fsmask |= FUSION_FLOWDEPTH;
                    }else{
                        fgr.fsmask &= (~FUSION_FLOWDEPTH);
                    }
                }
            }else if(fgr.fsmask & FUSION_FHOTGENE){// fusion not in public database
                if(srv){
                    if((srv < fuseOpt->mUsualFilter.mMinSupport && dpv < fuseOpt->mUsualFilter.mMinSupport) || 
                       (srs < fuseOpt->mUsualFilter.mMinSRSeed && dps < fuseOpt->mUsualFilter.mMinDPSeed)){
                        fgr.fsmask |= FUSION_FLOWSUPPORT;
                    }else{
                        fgr.fsmask &= (~FUSION_FLOWSUPPORT);
                    }
                }else{
                    fgr.fsmask |= FUSION_FLOWSUPPORT;
                }
                if(af < fuseOpt->mUsualFilter.mMinVAF){
                    fgr.fsmask |= FUSION_FLOWAF;
                }else{
                    fgr.fsmask &= (~FUSION_FLOWAF);
                }
                if(srv){
                    if((srv + srr) < fuseOpt->mUsualFilter.mMinDepth){
                        fgr.fsmask |= FUSION_FLOWDEPTH;
                    }else{
                        fgr.fsmask &= (~FUSION_FLOWDEPTH);
                    }
                }else{
                    fgr.fsmask |= FUSION_FLOWDEPTH;
                }
            }
            fgr.fusegene = fgl[i].hgene + "->" + fgl[i].tgene;                   // FusionGene
            if(fgl[i].hfrom1) fgr.fusepattern.append(trsl1[fgl[i].hidx].strand); // FusionPattern
            else fgr.fusepattern.append(trsl2[fgl[i].hidx].strand);
            if(fgl[i].tfrom1) fgr.fusepattern.append(trsl1[fgl[i].tidx].strand);
            else fgr.fusepattern.append(trsl2[fgl[i].tidx].strand);
            float sraf = 0, dpaf = 0;
            if(srr + srv) sraf = (double)(srv)/(double)(srv + srr);
            if(dpr + dpv) dpaf = (double)(dpv)/(double)(dpv + dpr);
            if(sraf > dpaf){                                                     // FusionReads TotalReads
                fgr.fusionreads = srv;
                fgr.totalreads = srr + srv;
            }else{
                fgr.fusionreads = dpv;
                fgr.totalreads = dpr + dpv;
            }
            fgr.fuserate = af;                                                   // FusionRate
            if(fgl[i].hfrom1){
                fgr.jctpos1 = start;                                             // JunctionPosition1
                fgr.strand1 = trsl1[fgl[i].hidx].strand;                         // Strand1;
                fgr.transcript1 = trsl1[fgl[i].hidx].getTrsWithVer();            // Transcript2
                fgr.exon1 = trsl1[fgl[i].hidx].exon;                             // exon1
            }else{
                fgr.jctpos1 = end;                                               // JunctionPosition1
                fgr.strand1 = trsl2[fgl[i].hidx].strand;                         // Strand1;
                fgr.transcript1 = trsl2[fgl[i].hidx].getTrsWithVer();            // Transcript1
                fgr.exon1 = trsl2[fgl[i].hidx].exon;                             // exon1
            }
            if(fgl[i].tfrom1){
                fgr.jctpos2 = start;                                             // JunctionPosition2
                fgr.strand2 = trsl1[fgl[i].tidx].strand;                         // Strand2
                fgr.transcript2 = trsl1[fgl[i].tidx].getTrsWithVer();            // Transcript2
                fgr.exon2 = trsl1[fgl[i].tidx].exon;                             // exon2
            }else{
                fgr.jctpos2 = end;                                               // JunctionPosition2
                fgr.strand2 = trsl2[fgl[i].tidx].strand;                         // Strand2
                fgr.transcript2 = trsl2[fgl[i].tidx].getTrsWithVer();            // Transcript2
                fgr.exon2 = trsl2[fgl[i].tidx].exon;                             // exon2
            }
            fgr.fusionsequence = vstr[16];                                       // FusionSequence
            fgr.fseqbp = std::atoi(vstr[17].c_str());                            // fseqBp
            fgr.indb = ((fgr.fsmask & FUSION_FINDB) ? "Y" : "N");                // inDB
            fgr.svt = vstr[0];                                                   // svType
            fgr.svsize = std::atoi(vstr[1].c_str());                             // svSize
            fgr.srcount = std::atoi(vstr[7].c_str());                            // srCount
            fgr.dpcount = std::atoi(vstr[8].c_str());                            // dpCount
            fgr.srrescued = srv;                                                 // srRescued
            fgr.dprescued = std::atoi(vstr[10].c_str());                         // dpRescued
            fgr.srrefcount  = std::atoi(vstr[11].c_str());                       // srRefCount
            fgr.dprefcount = std::atoi(vstr[12].c_str());                        // dpRefCount
            fgr.insbp = std::atoi(vstr[14].c_str());                             // insBp
            fgr.insseq = vstr[15];                                               // insSeq
            fgr.svid = std::atoi(vstr[18].c_str());                              // svID
            fgr.svint = std::atoi(vstr[19].c_str());                             // svInt
            fgr.fsHits = std::atoi(vstr[24].c_str());                            // fsHits;
            fgr.distance = fuseOpt->geneNear(fgr.gene1, fgr.chr1, fgr.jctpos1, fgr.gene2, fgr.chr2); // fpDist
            // mask FUSION_FINREPORTRNG now
            if(fuseOpt->mFsRptList.empty()){
                fgr.fsmask |= FUSION_FINREPORTRNG;
            }else{
                if(fuseOpt->inFsRptRange(fgr.gene1, fgr.gene2, fgr.exon1, fgr.exon2, "ee")){
                    fgr.fsmask |= FUSION_FINREPORTRNG;
                }else{
                    fgr.fsmask &= (~FUSION_FINREPORTRNG);
                }
            }
            if(fgr.fsmask & (FUSION_FINDB | FUSION_FMIRRORINDB)){
                if(fgr.fsHits >= 0 && fgr.fsHits <= fuseOpt->mWhiteFilter.mMaxRepHit) fgr.fsmask |= FUSION_FREALNPASSED;
                else fgr.fsmask &= (~FUSION_FREALNPASSED);
            }else{
                if(fgr.fsHits >= 0 && fgr.fsHits <= fuseOpt->mUsualFilter.mMaxRepHit) fgr.fsmask |= FUSION_FREALNPASSED;
                  else fgr.fsmask &= (~FUSION_FREALNPASSED);
            }
            if(fgr.fsmask & FUSION_DROP_MASK){
                fgr.fsmask &= (~(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
            }else{
                if((fgr.fsmask & FUSION_KEEP_MASK) != FUSION_KEEP_MASK){
                    fgr.fsmask &= (~(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
                }else{
                    if((fgr.fsmask & PRIMARY_KEEP_MASK) == PRIMARY_KEEP_MASK){
                        fgr.fsmask |= FUSION_FPRIMARY;
                    }else{
                        fgr.fsmask |= FUSION_FSUPPLEMENTARY;
                    }
                }
            }
            if((fgr.fsmask & fuseOpt->mFsMaskInclude) != fuseOpt->mFsMaskInclude){
                fgr.fsmask &= (~(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
            }
            if(fgr.fsmask & fuseOpt->mFsMaskExclude){
                fgr.fsmask &= (~(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
            }
            if(fgr.fsmask & FUSION_FCALLFROMRNASEQ){
                fgr.ts1name = vstr[25];                                            // ts1Name
                fgr.ts1pos = std::atoi(vstr[26].c_str());                          // ts1Pos
                fgr.ts2name = vstr[27];                                            // ts2Name
                fgr.ts2pos = std::atoi(vstr[28].c_str());                          // ts2Pos
            }
            if(debug){
                std::cout << fgr << std::endl;
            }
            frl.push_back(fgr);
        }
        bool reported = false;
        for(uint32_t fi = 0; fi < frl.size(); ++fi){
            if(frl[fi].fsmask & FUSION_FPRIMARY){
                fsrl.push_back(frl[fi]);
                reported = true;
            }
        }
        if(!reported){
            for(uint32_t fi = 0; fi < frl.size(); ++fi){
                if(frl[fi].fsmask & FUSION_FSUPPLEMENTARY){
                    fsrl.push_back(frl[fi]);
                }
            }
        }
        if(!fuseOpt->mSVModFile.empty()){
            std::vector<std::string> nfmsk;
            for(uint32_t fk = 0; fk < frl.size(); ++fk){
                nfmsk.push_back(std::to_string(frl[fk].fsmask));
            }
            vstr[23] = util::join(nfmsk, ";");
            std::string newrec;
            util::join(vstr, newrec, "\t");
            fsv << newrec << "\n";
        }
    }
    if(!fuseOpt->mSVModFile.empty()) fsv.close();
    fr.close();
}
