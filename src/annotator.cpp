#include "annotator.h"
#include "ThreadPool.h"

std::set<std::pair<int32_t, int32_t>>::iterator Annotator::getFirstOverlap(const std::set<std::pair<int32_t, int32_t>>& s, const std::pair<int32_t, int32_t>&p){
    auto itsv= std::upper_bound(s.begin(), s.end(), p);
    auto upIter = itsv;
    if(itsv != s.begin()){
        auto lowIter = --itsv;
        bool lowGet = false;
        while(p.first <= lowIter->second){
            if(itsv == s.begin()){
                return itsv;
            }
            --lowIter;
            lowGet = true;
        }
        if(lowGet){
            return ++lowIter;
        }
    }
    bool upGet = false;
    while(upIter != s.end() && p.second >= upIter->first){
        ++upIter;
        upGet = true;
    }
    if(upGet){
        return --upIter;
    }
    return s.end();
}

void Annotator::cgrsplit(const cgranges_t* cr, std::vector<RegItemCnt>& ctgRng, int32_t us){
    ctgRng.clear();
    int32_t ltid = -1, lbeg = -1, lend = -1, lsum = 0;
    for(int64_t i = 0; i < cr->n_r; ++i){
        const cr_intv_t *r = &cr->r[i];
        const char* name = cr->ctg[r->x>>32].name;
        int32_t tid = bam_name2id(mOpt->bamheader, name);
        if(r->label > us){// jut split these regions independently
            int32_t nsp = r->label / us + 1;
            std::vector<std::pair<int32_t, int32_t>> vpidx;
            util::divideVecIdx((int32_t)r->y - (int32_t)r->x, nsp, vpidx);
            for(uint32_t ki = 0; ki < vpidx.size(); ++ki){
                RegItemCnt irc;
                irc.mBeg = (int32_t)r->x + vpidx[ki].first;
                irc.mEnd = (int32_t)r->x + vpidx[ki].second;
                irc.mTid = tid;
                irc.mCount = us;
                if(ki) irc.mInterleved = true;
                ctgRng.push_back(irc);
            }
            if(lsum > 0){// append previous left ones
                RegItemCnt crc;
                crc.mBeg = lbeg;
                crc.mEnd = lend;
                crc.mTid = ltid;
                crc.mCount = lsum;
                ctgRng.push_back(crc);
            }
            ltid = lbeg = lend = -1; lsum = 0; // reinitialize
            continue;
        }
        if(tid == ltid){
            if(lsum >= us){
                RegItemCnt crc;
                crc.mBeg = lbeg;
                crc.mEnd = lend;
                crc.mTid = ltid;
                crc.mCount = lsum;
                ctgRng.push_back(crc);
                lbeg = (int32_t)r->x;
                lend = (int32_t)r->y;
                lsum = r->label;
            }else{
                lsum += r->label;
                lend = (int32_t)r->y;
            }
        }else{
            if(ltid > 0){
                RegItemCnt crc;
                crc.mBeg = lbeg;
                crc.mEnd = lend;
                crc.mTid = ltid;
                crc.mCount = lsum;
                ctgRng.push_back(crc);
            }
            ltid = tid;
            lbeg = (int32_t)r->x;
            lend = (int32_t)r->y;
            lsum = r->label;
        }
    }
    if(lsum > 0){
        RegItemCnt crc;
        crc.mTid = ltid;
        crc.mBeg = lbeg;
        crc.mEnd = lend;
        crc.mCount = lsum;
        ctgRng.push_back(crc);
    }
}

Stats* Annotator::covAnnotate(std::vector<SVRecord>& svs){
    // Store all regions of SV into cgranges_t
    util::loginfo("Beg construct SVs cgranges_t");
    cgranges_t* crsv = cr_init();
    int32_t regBeg = -1, regEnd = -1;
    for(uint32_t i = 0; i < svs.size(); ++i){
        regBeg = std::max(0, svs[i].mSVStart - mOpt->libInfo->mMaxNormalISize);
        regEnd = std::min(svs[i].mSVStart + mOpt->libInfo->mMaxNormalISize, (int32_t)mOpt->bamheader->target_len[svs[i].mChr1]);
        cr_add(crsv, svs[i].mNameChr1.c_str(), regBeg, regEnd, 1);
        regBeg = std::max(0, svs[i].mSVEnd - mOpt->libInfo->mMaxNormalISize);
        regEnd = std::min(svs[i].mSVEnd + mOpt->libInfo->mMaxNormalISize, (int32_t)mOpt->bamheader->target_len[svs[i].mChr2]);
        cr_add(crsv, svs[i].mNameChr2.c_str(), regBeg, regEnd, 1);
    }
    if(!cr_is_sorted(crsv)) cr_sort(crsv);
    cr_merge_pre_index(crsv);
    std::vector<cgranges_t*> svregs(mOpt->contigNum, NULL);
    for(uint32_t i = 0; i < svregs.size(); ++i) svregs[i] = cr_init();
    for(int64_t i = 0; i < crsv->n_r; ++i){
        const cr_intv_t *r = &crsv->r[i];
        const char* name = crsv->ctg[r->x>>32].name;
        int32_t tid = bam_name2id(mOpt->bamheader, name);
        cr_add(svregs[tid], name, (int32_t)r->x, (int32_t)r->y, r->label);
    }
    for(uint32_t i = 0; i < svregs.size(); ++i) cr_index2(svregs[i], 1);
    util::loginfo("End construct SVs cgranges_t");
    util::loginfo("Beg split SVs cgranges_t");
    std::vector<RegItemCnt> ctgRng;
    cgrsplit(crsv, ctgRng, 1000);
    std::sort(ctgRng.begin(), ctgRng.end());
    cr_destroy(crsv);
    util::loginfo("End split Svs cgranges_t, got: " + std::to_string(ctgRng.size()) + " sub regions");
    // Preprocess REF and ALT
    ContigBpRegions bpRegion(mOpt->bamheader->n_targets);
    std::vector<int32_t> refAlignedReadCount(svs.size());
    std::vector<int32_t> refAlignedSpanCount(svs.size());
    util::loginfo("Beg extracting breakpoint regions of each precisely classified SV");
    for(auto itsv = svs.begin(); itsv != svs.end(); ++ itsv){
        if(!itsv->mPrecise) continue;
        // Iterate all break point
        for(int bpPoint = 0; bpPoint < 2; ++bpPoint){
            int32_t regChr, regStart, regEnd, bpPos;
            if(bpPoint){// SV ending position region
                regChr = itsv->mChr2;
                regStart = std::max(0, itsv->mSVEnd - mOpt->filterOpt->mMinFlankSize);
                regEnd = std::min(itsv->mSVEnd + mOpt->filterOpt->mMinFlankSize, (int32_t)mOpt->bamheader->target_len[itsv->mChr2]);
                bpPos = itsv->mSVEnd;
            }else{// SV starting position region
                regChr = itsv->mChr1;
                regStart = std::max(0, itsv->mSVStart - mOpt->filterOpt->mMinFlankSize);
                regEnd = std::min(itsv->mSVStart + mOpt->filterOpt->mMinFlankSize, (int32_t)mOpt->bamheader->target_len[itsv->mChr1]);
                bpPos = itsv->mSVStart;
            }
            BpRegion br;
            br.mIsSVEnd = bpPoint;
            br.mBpPos = bpPos;
            br.mID = itsv->mID;
            br.mRegEnd = regEnd;
            br.mRegStart = regStart;
            br.mSVT = itsv->mSVT;
            bpRegion[regChr].push_back(br);
        }
    }
    // Sort all BpRegion by breakpoint position on reference
    for(auto& refIndex : mOpt->svRefID) std::sort(bpRegion[refIndex].begin(), bpRegion[refIndex].end());
    util::loginfo("End extracting breakpoint regions of each precisely classified SV");
    // Get spanning point regions
    util::loginfo("Beg extracting PE supported breakpoints of each SV");
    ContigSpanPoints spanPoint;
    spanPoint.resize(mOpt->contigNum);
    for(auto itsv = svs.begin(); itsv != svs.end(); ++itsv){
        if(itsv->mPESupport == 0) continue;
        spanPoint[itsv->mChr1].push_back(SpanPoint(itsv->mSVStart, itsv->mSVT, itsv->mID, false));
        spanPoint[itsv->mChr2].push_back(SpanPoint(itsv->mSVEnd, itsv->mSVT, itsv->mID, true));
    }
    for(uint32_t i = 0; i < spanPoint.size(); ++i) std::sort(spanPoint[i].begin(), spanPoint[i].end());
    util::loginfo("End extracting PE supported breakpoints of each SV");
    // Get coverage from each contig in parallel
    Stats* covStats = new Stats(mOpt, svs.size());
    std::vector<std::future<void>> statRets(ctgRng.size());
    for(uint32_t i = 0; i < ctgRng.size(); ++i){
        statRets[i] = mOpt->pool->enqueue(&Stats::stat, covStats, std::ref(svs),  std::ref(bpRegion), 
                                          std::ref(spanPoint), std::ref(ctgRng[i]), svregs[ctgRng[i].mTid]);
    }
    for(auto& e: statRets) e.get();
    for(auto& e: svregs) cr_destroy(e);
    return covStats;
}

void Annotator::getDNABpTrs(TrsRecList& trl, const std::string& chr, int32_t pos, htsFile* fp, tbx_t* tbx){
    std::vector<TrsRec> trsList;
    std::stringstream regs;
    regs << chr << ":" << pos << "-" << pos + 1;
    hts_itr_t* itr = tbx_itr_querys(tbx, regs.str().c_str());
    kstring_t rec = {0, 0, 0};
    std::vector<std::string> vstr;
    std::set<std::string> trsGotIE; // gene got intron/exon spanning this bp
    std::map<std::string, std::string> gene2cnc; // gene with canonical transcript spanning this bp
    std::map<std::string, std::string> gene2rdt; // gene with no canonical transcript spanning this bp, keep random one
    while(tbx_itr_next(fp, tbx, itr, &rec) >= 0){
        util::split(rec.s, vstr, "\t");
        TrsRec tr;
        tr.strand = vstr[3];
        tr.unit = vstr[4];
        tr.number = vstr[5];
        tr.name = vstr[6];
        tr.gene = vstr[7];
        tr.version = vstr[8];
        tr.primary = vstr[9];
        tr.chr = vstr[0];
        tr.drop = false;
        tr.pos = pos;
        trsList.push_back(tr);
        if(!util::startsWith(tr.unit, "utr")) trsGotIE.insert(tr.name);
        if(tr.primary == "Y") gene2cnc[tr.gene] = tr.name;
    }
    // get clean transcript list spanning this breakpoint
    for(auto& tr: trsList){
        // if utr and exon/intron of one trnascript spanning this breakpoint simultaneously, keep exon and intron only
        if(trsGotIE.find(tr.name) != trsGotIE.end()){
            if(util::startsWith(tr.unit, "utr")){
                tr.drop = true;
            }
        }
        // if one gene has canonical transcript spanning this breakpoint, keep the canonical transcript only
        // if one gene has none canonical transcript spanning this breakpoint, random choose one trascript of this gene
        auto g2citr = gene2cnc.find(tr.gene);
        if(g2citr != gene2cnc.end()){
            if(tr.name != g2citr->second) tr.drop = true;
        }else{
            if(!tr.drop) gene2rdt[tr.gene] = tr.name;
        }
    }
    for(auto& tr: trsList){
        if(tr.drop) continue;
        auto g2rtr = gene2rdt.find(tr.gene);
        if(g2rtr != gene2rdt.end()){
            if(tr.name != g2rtr->second){
                tr.drop = true;
            }
        }
    }
    // output valid transcript
    for(auto& tr: trsList){
        if(!tr.drop) trl.push_back(tr);
    }
    // cleanup
    if(itr) hts_itr_destroy(itr);
    if(rec.s) free(rec.s);
}

void Annotator::geneAnnoDNA(SVSet& svs, GeneInfoList& gl){
    gl.resize(svs.size());
    // split range list 
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    util::divideVecIdx(svs.size(), mOpt->nthread, vpidx);
    // parallel run
    std::vector<std::future<void>> annRets(vpidx.size());
    for(uint32_t i = 0; i < vpidx.size(); ++i){
        annRets[i] = mOpt->pool->enqueue(&Annotator::rangeGeneAnnoDNA, this, std::ref(svs), std::ref(gl), vpidx[i].first, vpidx[i].second);
    }
    for(auto& e: annRets) e.get();
}

void Annotator::rangeGeneAnnoDNA(SVSet& svs, GeneInfoList& gl, int32_t begIdx, int32_t endIdx){
    htsFile* fp = hts_open(mOpt->annodb.c_str(), "r");
    tbx_t* tbx = tbx_index_load(mOpt->annodb.c_str());
    for(int32_t i = begIdx; i < endIdx; ++i){
        // get trascripts at breakpoint 1
        getDNABpTrs(gl[i].mGene1, svs[i].mNameChr1, svs[i].mSVStart, fp, tbx);
        gl[i].mPos1 = svs[i].mSVStart;
        gl[i].mChr1 = svs[i].mNameChr1;
        // get transcripts at breakpoint 2
        getDNABpTrs(gl[i].mGene2, svs[i].mNameChr2, svs[i].mSVEnd, fp, tbx);
        gl[i].mPos2 = svs[i].mSVEnd;
        gl[i].mChr2 = svs[i].mNameChr2;
        // annotate fusion gene
        for(uint32_t g1 = 0; g1 < gl[i].mGene1.size(); ++g1){
            for(uint32_t g2 = 0; g2 < gl[i].mGene2.size(); ++g2){
                FuseGene fsg = svutil::getFusionGene(gl[i].mGene1[g1].gene, gl[i].mGene2[g2].gene, gl[i].mGene1[g1].strand[0], gl[i].mGene2[g2].strand[0], svs[i].mSVT);
                if(fsg.status & FUSION_FHTFLSWAPPED){
                    // remember h/t gene sources
                    fsg.hfrom1 = false;
                    fsg.hidx = g2;
                    fsg.tfrom1 = true;
                    fsg.tidx = g1;
                    // add cigar string of catentaion around breakpoint gl[i].mGene2[g2] -> gl[i].mGene1[g1], TODO...
                }else{
                    // remember h/t gene sources
                    fsg.hfrom1 = true;
                    fsg.hidx = g1;
                    fsg.tfrom1 = false;
                    fsg.tidx = g2;
                    // add cigar string of catentaion around breakpoint gl[i].mGene1[g1] -> gl[i].mGene2[g2], TODO...
                }
                gl[i].mFuseGene.push_back(fsg);
            }
        }
    }
    if(tbx) tbx_destroy(tbx);
    hts_close(fp);
}

void Annotator::getRNABpTrs(TrsRecList& trl, const std::string& chr, int32_t pos, htsFile* fp, tbx_t* tbx){
    std::stringstream regs;
    regs << chr << ":" << pos << "-" << pos + 1;
    hts_itr_t* itr = tbx_itr_querys(tbx, regs.str().c_str());
    kstring_t rec = {0, 0, 0};
    std::vector<std::string> vstr;
    std::set<std::string> trsGotIE; // gene got intron/exon spanning this bp
    TrsRecList trsList;
    while(tbx_itr_next(fp, tbx, itr, &rec) >= 0){
        util::split(rec.s, vstr, "\t");
        TrsRec tr;
        tr.strand = vstr[9];
        tr.version = vstr[10];
        tr.unit = vstr[3];
        tr.number = vstr[4];
        tr.name = vstr[0];
        tr.gene = vstr[5];
        tr.chr = vstr[6];
        tr.drop = false;
        tr.pos = svutil::trpos2gnpos(pos, std::atoi(vstr[1].c_str()), std::atoi(vstr[2].c_str()), std::atoi(vstr[7].c_str()),  vstr[9][0]);
        if(tr.strand[0] == '+'){
            tr.eoffset = std::atoi(vstr[8].c_str()) - tr.pos;
            tr.ioffset = tr.pos - std::atoi(vstr[7].c_str());
        }else{
            tr.eoffset = tr.pos - std::atoi(vstr[7].c_str());
            tr.ioffset = std::atoi(vstr[8].c_str()) - tr.pos;
        }
        trsList.push_back(tr);
        if(!util::startsWith(tr.unit, "utr")) trsGotIE.insert(tr.name);
    }
    // get clean transcript list spanning this breakpoint
    for(auto& tr: trsList){
        // if utr and exon/intron of one trnascript spanning this breakpoint simultaneously, keep exon and intron only
        if(trsGotIE.find(tr.name) != trsGotIE.end()){
            if(util::startsWith(tr.unit, "utr")){
                tr.drop = true;
            }
        }
    }
    // output valid transcript
    for(auto& tr: trsList){
        if(!tr.drop) trl.push_back(tr);
    }
    // cleanup
    if(itr) hts_itr_destroy(itr);
    if(rec.s) free(rec.s);
}

void Annotator::geneAnnoRNA(SVSet& svs, GeneInfoList& gl){
    gl.resize(svs.size());
    // split range list 
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    util::divideVecIdx(svs.size(), mOpt->nthread, vpidx);
    // parallel run
    std::vector<std::future<void>> annRets(vpidx.size());
    for(uint32_t i = 0; i < vpidx.size(); ++i){
        annRets[i] = mOpt->pool->enqueue(&Annotator::rangeGeneAnnoRNA, this, std::ref(svs), std::ref(gl), vpidx[i].first, vpidx[i].second);
    }
    for(auto& e: annRets) e.get();
}

void Annotator::rangeGeneAnnoRNA(SVSet& svs, GeneInfoList& gl, int32_t begIdx, int32_t endIdx){
    htsFile* fp = hts_open(mOpt->annodb.c_str(), "r");
    tbx_t* tbx = tbx_index_load(mOpt->annodb.c_str());
    for(int32_t i = begIdx; i < endIdx; ++i){
        // get trascript at breakpoint 1
        getRNABpTrs(gl[i].mGene1, svs[i].mNameChr1, svs[i].mSVStart, fp, tbx);
        gl[i].mChr1 = gl[i].mGene1[0].chr;
        gl[i].mPos1 = gl[i].mGene1[0].pos;
        // get trascript at breakpoint 2
        getRNABpTrs(gl[i].mGene2, svs[i].mNameChr2, svs[i].mSVEnd, fp, tbx);
        gl[i].mChr2 = gl[i].mGene2[0].chr;
        gl[i].mPos2 = gl[i].mGene2[0].pos;
        // annotate fusion gene
        for(uint32_t g1 = 0; g1 < gl[i].mGene1.size(); ++g1){
            for(uint32_t g2 = 0; g2 < gl[i].mGene2.size(); ++g2){
                FuseGene fsg = svutil::getFusionGene(gl[i].mGene1[g1].gene, gl[i].mGene2[g2].gene, '+', '+', svs[i].mSVT);
                if(fsg.status & FUSION_FHTFLSWAPPED){
                    // remember h/t gene sources
                    fsg.hfrom1 = false;
                    fsg.hidx = g2;
                    fsg.tfrom1 = true;
                    fsg.tidx = g1;
                    // add cigar string of catentaion around breakpoint gl[i].mGene2[g2] -> gl[i].mGene1[g1]
                    fsg.cigar = svutil::bp2cigar(gl[i].mGene2[g2], gl[i].mGene1[g1], svs[i].mSVT);
                }else{
                    // remember h/t gene sources
                    fsg.hfrom1 = true;
                    fsg.hidx = g1;
                    fsg.tfrom1 = false;
                    fsg.tidx = g2;
                    // add cigar string of catentaion around breakpoint gl[i].mGene1[g1] -> gl[i].mGene2[g2]
                    fsg.cigar = svutil::bp2cigar(gl[i].mGene1[g1], gl[i].mGene2[g2], svs[i].mSVT);
                }
                gl[i].mFuseGene.push_back(fsg);
            }
        }
    }
    if(tbx) tbx_destroy(tbx);
    hts_close(fp);
}
