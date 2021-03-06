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

Stats* Annotator::covAnnotate(SVSet& svs){
    // Store all regions of SV into cgranges_t
    util::loginfo("Beg construct SVs cgranges_t");
    cgranges_t* crsv = cr_init();
    int32_t regBeg = -1, regEnd = -1;
    for(uint32_t i = 0; i < svs.size(); ++i){
        regBeg = std::max(0, svs[i]->mSVStart - mOpt->libInfo->mMaxNormalISize);
        regEnd = std::min(svs[i]->mSVStart + mOpt->libInfo->mMaxNormalISize, (int32_t)mOpt->bamheader->target_len[svs[i]->mChr1]);
        cr_add(crsv, svs[i]->mNameChr1.c_str(), regBeg, regEnd, 1);
        regBeg = std::max(0, svs[i]->mSVEnd - mOpt->libInfo->mMaxNormalISize);
        regEnd = std::min(svs[i]->mSVEnd + mOpt->libInfo->mMaxNormalISize, (int32_t)mOpt->bamheader->target_len[svs[i]->mChr2]);
        cr_add(crsv, svs[i]->mNameChr2.c_str(), regBeg, regEnd, 1);
    }
    if(!cr_is_sorted(crsv)) cr_sort(crsv);
    cr_merge_pre_index(crsv);
#ifdef DEBUG
    if(mOpt->debug & DEBUG_FANNC){
        std::cout << "DEBUG_ANNO_CGRANGE_OF_SVS: " << std::endl;
        cr_iter_usual(crsv, stdout);
    }
#endif
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
#ifdef DEBUG
    if(mOpt->debug & DEBUG_FANNC){
        std::cout << "DEBUG_SV_CGRANGE_CONSTRUCTED: " << std::endl;
        for(uint32_t i = 0; i < svregs.size(); ++i){
            cr_iter_indexed(svregs[i], stdout);
        }
    }
#endif
    util::loginfo("Beg split SVs cgranges_t");
    std::vector<RegItemCnt> ctgRng;
    cgrsplit(crsv, ctgRng, mOpt->batchsvn);
    std::sort(ctgRng.begin(), ctgRng.end());
#ifdef DEBUG
    if(mOpt->debug & DEBUG_FANNC){
        std::cout << "DEBUG_SV_CONTIG_RANGE_LIST: " << std::endl;
        for(uint32_t i = 0; i < ctgRng.size(); ++i){
            std::cout << ctgRng[i] << std::endl;
        }
    }
#endif
    cr_destroy(crsv);
    util::loginfo("End split Svs cgranges_t, got: " + std::to_string(ctgRng.size()) + " sub regions");
    // Preprocess REF and ALT
    ContigBpRegions bpRegion(mOpt->bamheader->n_targets);
    std::vector<int32_t> refAlignedReadCount(svs.size());
    std::vector<int32_t> refAlignedSpanCount(svs.size());
    util::loginfo("Beg extracting breakpoint regions of each precisely classified SV");
    for(auto& itsv: svs){
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
    for(auto& itsv: svs){
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
    // if empty, append an empty tr
    if(trl.empty()){
        TrsRec tmpTrs;
        tmpTrs.chr = chr;
        tmpTrs.pos = pos;
        trl.push_back(tmpTrs);
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
        getDNABpTrs(gl[i].mGene1, svs[i]->mNameChr1, svs[i]->mSVStart, fp, tbx);
        gl[i].mPos1 = svs[i]->mSVStart;
        gl[i].mChr1 = svs[i]->mNameChr1;
        // get transcripts at breakpoint 2
        getDNABpTrs(gl[i].mGene2, svs[i]->mNameChr2, svs[i]->mSVEnd, fp, tbx);
        gl[i].mPos2 = svs[i]->mSVEnd;
        gl[i].mChr2 = svs[i]->mNameChr2;
        // annotate fusion gene
        std::set<std::string> fgAdded;
        for(uint32_t g1 = 0; g1 < gl[i].mGene1.size(); ++g1){
            for(uint32_t g2 = 0; g2 < gl[i].mGene2.size(); ++g2){
                svutil::getexon(gl[i].mGene1[g1], gl[i].mGene2[g2], svs[i]->mSVT);
                FuseGene fsg = svutil::getFusionGene(gl[i].mGene1[g1].gene, gl[i].mGene2[g2].gene, gl[i].mGene1[g1].strand[0], gl[i].mGene2[g2].strand[0], svs[i]->mSVT);
                gl[i].mGene1[g1].getCatPart(svs[i]->mSVT, true);
                gl[i].mGene2[g2].getCatPart(svs[i]->mSVT, false);
#ifdef DEBUG
                if(mOpt->debug & DEBUG_FANNG){
                    std::cout << gl[i] << std::endl;
                    std::cout << fsg.debugStr() << std::endl;
                }
#endif
                std::string newFg = fsg.hgene + "->" + fsg.tgene;
                if(fgAdded.find(newFg) != fgAdded.end()) continue; // do not add twice for A,B->A,B
                fgAdded.insert(newFg);
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

void Annotator::getRNABpTrs(TrsRecList& trl, const std::string& chr, int32_t pos, htsFile* fp, tbx_t* tbx, bool isbp1, int32_t svt){
    hts_itr_t* itr= tbx_itr_querys(tbx, chr.c_str());
    kstring_t rec = {0, 0, 0};
    std::vector<std::string> vstr;
    // first run get all exon units of this transcript
    std::vector<Rna2DnaUnit> r2dl;
    while(tbx_itr_next(fp, tbx, itr, &rec) >= 0){
        util::split(rec.s, vstr, "\t");
        if(util::startsWith(vstr[3], "utr")) continue;
        Rna2DnaUnit r2du;
        r2du.tname = vstr[0];
        r2du.tbeg = std::atoi(vstr[1].c_str());
        r2du.tend = std::atoi(vstr[2].c_str());
        r2du.uname = vstr[3];
        r2du.ucount = std::atoi(vstr[4].c_str());
        r2du.gname = vstr[5];
        r2du.gchr = vstr[6];
        r2du.gbeg = std::atoi(vstr[7].c_str());
        r2du.gend = std::atoi(vstr[8].c_str());
        r2du.gstrand = vstr[9][0];
        r2du.tversion = vstr[10];
        r2dl.push_back(r2du);
    }
    std::sort(r2dl.begin(), r2dl.end());
    TrsRec tr;
    svutil::getPropTrs(r2dl, pos, svt, isbp1, tr);
    trl.push_back(tr);
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
        getRNABpTrs(gl[i].mGene1, svs[i]->mNameChr1, svs[i]->mSVStart, fp, tbx, true, svs[i]->mSVT);
        gl[i].mChr1 = gl[i].mGene1[0].chr;
        gl[i].mPos1 = gl[i].mGene1[0].pos;
        // get trascript at breakpoint 2
        getRNABpTrs(gl[i].mGene2, svs[i]->mNameChr2, svs[i]->mSVEnd, fp, tbx, false, svs[i]->mSVT);
        gl[i].mChr2 = gl[i].mGene2[0].chr;
        gl[i].mPos2 = gl[i].mGene2[0].pos;
        // annotate fusion gene
        for(uint32_t g1 = 0; g1 < gl[i].mGene1.size(); ++g1){
            for(uint32_t g2 = 0; g2 < gl[i].mGene2.size(); ++g2){
                FuseGene fsg = svutil::getFusionGene(gl[i].mGene1[g1].gene, gl[i].mGene2[g2].gene, '+', '+', svs[i]->mSVT);
                gl[i].mGene1[g1].getCatPart(svs[i]->mSVT, true);
                gl[i].mGene2[g2].getCatPart(svs[i]->mSVT, false);
                gl[i].mGene1[g1].getCigar();
                gl[i].mGene2[g2].getCigar();
                if(fsg.status & FUSION_FHTFLSWAPPED){
                    // remember h/t gene sources
                    fsg.hfrom1 = false;
                    fsg.hidx = g2;
                    fsg.tfrom1 = true;
                    fsg.tidx = g1;
                }else{
                    // remember h/t gene sources
                    fsg.hfrom1 = true;
                    fsg.hidx = g1;
                    fsg.tfrom1 = false;
                    fsg.tidx = g2;
                }
                gl[i].mFuseGene.push_back(fsg);
            }
        }
    }
    if(tbx) tbx_destroy(tbx);
    hts_close(fp);
}

void Annotator::refineCovAnno(Stats* sts, const SVSet& svs){
    if(mOpt->bamout.empty()) return;
    ReadSupportStatMap rssm;
    PePtnMap pem;
    util::loginfo("Beg collect SV supporting reads info");
    getReadSupportStatus(mOpt->bamout, svs, rssm, pem, mOpt->realnf);
    util::loginfo("End collect SV supporting reads info");
    util::loginfo("Beg collect PE partner reads");
    // temporary out bam
    std::string tmpBam = mOpt->bamout + ".tmp.bam";
    samFile *tmpSamFp = sam_open(tmpBam.c_str(), "wb");
    // rescue and stat PE
    BedRegs *br = new BedRegs();
    br->mCR = cr_init();
    for(auto& e: pem){
        cr_add(br->mCR, 
               e.second->chr.c_str(),
               std::max(0, e.second->mpos - mOpt->libInfo->mReadLen), 
               e.second->mpos + mOpt->libInfo->mReadLen, 0);
    }
    cr_index2(br->mCR, 1);
    samFile *ttSamFp = sam_open(mOpt->bamfile.c_str(), "r");
    bam_hdr_t *h = sam_hdr_read(ttSamFp);
    hts_idx_t *idx = sam_index_load(ttSamFp, mOpt->bamfile.c_str());
    cgranges_t *cr = br->mCR;
    assert(sam_hdr_write(tmpSamFp, h) >= 0);
    bam1_t *b = bam_init1();
    for(int32_t ctg_id = 0; ctg_id < cr->n_ctg; ++ctg_id){
        int64_t i, *xb = 0, max_b = 0, n = 0;
        n = cr_overlap_int(cr, ctg_id, 0, INT_MAX, &xb, &max_b);
        for(i = 0; i < n; ++i){
            const char* chr = cr->ctg[ctg_id].name;
            int32_t beg = cr_start(cr, xb[i]), end = cr_end(cr, xb[i]);
            hts_itr_t *itr = sam_itr_queryi(idx, sam_hdr_name2tid(h, chr), beg, end);
            while(sam_itr_next(ttSamFp, itr, b) >= 0){
                if(b->core.tid < b->core.mtid || (b->core.tid == b->core.mtid && b->core.pos <= b->core.mpos)){
                    if(bam_aux_get(b, "SA")) continue; // skip
                    auto iter = pem.find(bam_get_qname(b));
                    if(iter != pem.end()){
                        if((iter->second->is_read1 && b->core.flag & BAM_FREAD1) ||
                           (!iter->second->is_read1 && b->core.flag & BAM_FREAD2)){
                            if(!iter->second->found && b->core.qual > mOpt->filterOpt->minMapQual){
                                iter->second->valid = true;
                                iter->second->found = true;
                                bam_aux_update_int(b, "ZF", iter->second->svid);
                                bam_aux_update_int(b, "ST", 1);
                                assert(sam_write1(tmpSamFp, h, b) >= 0);
                            }
                        }
                    }
                }
            }
            sam_itr_destroy(itr);
        }
        free(xb);
    }
    hts_idx_destroy(idx);
    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(ttSamFp);
    delete br;
    // update mPEWithPtn
    for(auto& e: pem){
        if(e.second->valid) svs[e.second->svid]->mPEWithPtn += 1;
        delete e.second;
    }
    util::loginfo("End collect PE partner reads");
    util::loginfo("Beg resolve reads supporting multiple SV");
    std::map<std::string, int32_t> drec;
    // first run, resolve reads supporting multiple svs
    for(auto iter = rssm.begin(); iter != rssm.end(); ++iter){
        if(iter->second->mR1MapQ && iter->second->mR2MapQ){
            if(iter->second->mR1SVID != iter->second->mR2SVID){
                // find the one which is not in repeat region
                bool r1svrp = (svs[iter->second->mR1SVID]->mRealnRet < 0 || svs[iter->second->mR1SVID]->mRealnRet > mOpt->fuseOpt->mWhiteFilter.mMaxRepHit);
                bool r2svrp = (svs[iter->second->mR2SVID]->mRealnRet < 0 || svs[iter->second->mR2SVID]->mRealnRet > mOpt->fuseOpt->mWhiteFilter.mMaxRepHit);
                if(r1svrp == r2svrp){
                    if(!r1svrp){ // both not in repeat region
                        if(svs[iter->second->mR1SVID]->mSRSupport < svs[iter->second->mR2SVID]->mSRSupport){
                            r1svrp = true;
                            r2svrp = false;
                        }else if(svs[iter->second->mR1SVID]->mSRSupport == svs[iter->second->mR2SVID]->mSRSupport &&
                                 svs[iter->second->mR1SVID]->mPESupport < svs[iter->second->mR2SVID]->mPESupport){
                            r1svrp = true;
                            r2svrp = false;
                        }else{
                            r1svrp = false;
                            r2svrp = true;
                        }
                    }
                }
                if(r1svrp){
                    sts->mTotalAltCnts[iter->second->mR1SVID] -= 1;
                    if(iter->second->mR1SRT){
                        sts->mSpnCnts[iter->second->mR1SVID].mAltCnt -= 1;
                        sts->mSpnCnts[iter->second->mR1SVID].mAltQual[iter->second->mR1MapQ] -= 1;
                    }else{
                        sts->mJctCnts[iter->second->mR1SVID].mAltCnt -= 1;
                        sts->mJctCnts[iter->second->mR1SVID].mAltQual[iter->second->mR1MapQ] -= 1;
                    }
                }
                if(r2svrp){
                    sts->mTotalAltCnts[iter->second->mR2SVID] -= 1;
                    if(iter->second->mR2SRT){
                        sts->mSpnCnts[iter->second->mR2SVID].mAltCnt -= 1;
                        sts->mSpnCnts[iter->second->mR2SVID].mAltQual[iter->second->mR2MapQ] -= 1;
                    }else{
                        sts->mJctCnts[iter->second->mR2SVID].mAltCnt -= 1;
                        sts->mJctCnts[iter->second->mR2SVID].mAltQual[iter->second->mR2MapQ] -= 1;
                    }
                }
                if(r1svrp && r2svrp){
                    drec[iter->first] = 3;
                    iter->second->mR1Hit = 0;
                    iter->second->mR2Hit = 0;
                    iter->second->mR1Seed = 0;
                    iter->second->mR2Seed = 0;
                }else if(r1svrp){
                    drec[iter->first] = 1;
                    iter->second->mR1Hit = 0;
                    iter->second->mR1Seed = 0;
                }else if(r2svrp){
                    drec[iter->first] = 2;
                    iter->second->mR2Hit = 0;
                    iter->second->mR2Seed = 0;
                }
            }else{
                sts->mTotalAltCnts[iter->second->mR1SVID] -= 1;
            }
        }
    }
    util::loginfo("End resolve reads supporting multiple SV");
    util::loginfo("Beg estimate reads in repeat region");
    // second run, stat one end multiple mapping status
    if(mOpt->overlapRegs){
        std::map<int32_t, SeedStatus> mss;
        for(auto iter = rssm.begin(); iter != rssm.end(); ++iter){
            bool r1rpt = false , r2rpt = false;
            if(iter->second->mR1Hit > 2){
                // check r1
                if(iter->second->mR1PHit > 1){
                    if(!mOpt->overlapRegs->overlap(iter->second->mR1PChr.c_str(), iter->second->mR1PBeg, iter->second->mR1PEnd)){// not in probe region
                        r1rpt = true;
                        iter->second->mR1PTgt = 1;
                        iter->second->mR1Seed = 0;
                    }
                }
                if(!r1rpt){
                    if(iter->second->mR1SHit > 1){
                        if(!mOpt->overlapRegs->overlap(iter->second->mR1SChr.c_str(), iter->second->mR1SBeg, iter->second->mR1SEnd)){// not in probe region
                            r1rpt = true;
                            iter->second->mR1STgt = 1;
                            iter->second->mR1Seed = 0;
                        }
                    }
                }
            }
            if(iter->second->mR2Hit > 2){
                // check r2
                if(iter->second->mR2PHit > 1){
                    if(!mOpt->overlapRegs->overlap(iter->second->mR2PChr.c_str(), iter->second->mR2PBeg, iter->second->mR2PEnd)){// not in probe region
                        r2rpt = true;
                        iter->second->mR2PTgt = 1;
                        iter->second->mR2Seed = 0;
                    }
                }
                if(!r2rpt){
                    if(iter->second->mR2SHit > 1){
                        if(!mOpt->overlapRegs->overlap(iter->second->mR2SChr.c_str(), iter->second->mR2SBeg, iter->second->mR2SEnd)){// not in probe region
                            r2rpt = true;
                            iter->second->mR2STgt = 1;
                            iter->second->mR2Seed = 0;
                        }
                    }
                }
            }
            if(iter->second->mR1SVID >= 0){
                auto siter = mss.find(iter->second->mR1SVID);
                if(siter == mss.end()){
                    SeedStatus ss;
                    if(r1rpt){
                        ss.mmapsrt = 1;
                        ss.allsrt = 1;
                    }else{
                        ss.mmapsrt = 0;
                        ss.allsrt = iter->second->mR1Seed;
                    }
                    mss[iter->second->mR1SVID] = ss;
                }else{
                    if(r1rpt){
                        siter->second.mmapsrt += 1;
                        siter->second.allsrt += 1;
                    }else{
                        siter->second.allsrt += iter->second->mR1Seed;
                    }
                }
            }
            if(iter->second->mR2SVID >= 0){
                auto siter = mss.find(iter->second->mR2SVID);
                if(siter == mss.end()){
                    SeedStatus ss;
                    if(r2rpt){
                        ss.mmapsrt = 1;
                        ss.allsrt = 1;
                    }else{
                        ss.mmapsrt = 0;
                        ss.allsrt = iter->second->mR2Seed;
                    }
                    mss[iter->second->mR2SVID] = ss;
                }else{
                    if(r2rpt){
                        siter->second.mmapsrt += 1;
                        siter->second.allsrt += 1;
                    }else{
                        siter->second.allsrt += iter->second->mR2Seed;
                    }
                }
            }
        }
        // stat sr event rescued sr seed multiple mapping rate
        for(auto iter = mss.begin(); iter != mss.end(); ++iter){
            svs[iter->first]->mSRSResMAlnCnt = iter->second.mmapsrt;
            svs[iter->first]->mSRSResAllCnt = iter->second.allsrt;
        }
    }
    // third run, stat fusion read pattern
    for(auto iter = rssm.begin(); iter != rssm.end(); ++iter){
        int r1pt = -1, r2pt = -1;
        int svid = -1;
        if(iter->second->mR1SVID >= 0) svid = iter->second->mR1SVID;
        else if(iter->second->mR2SVID >= 0) svid =iter->second->mR2SVID;
        if(svid >= 0){
            iter->second->countPattern(svs[svid]->mChr1, svs[svid]->mSVStart, r1pt, r2pt, svs[svid]->mSRSupport > 0);
            if(r1pt >= 0) svs[svid]->mFsPattern[r1pt] += 1;
            if(r2pt >= 0) svs[svid]->mFsPattern[r2pt] += 1;
        }
    }
    util::loginfo("End estimate reads in repeat region");
    if(mOpt->refinedp){
        // forth run, update dprescue and dpcount
        for(uint32_t i = 0; i < svs.size(); ++i){
            if(svs[i]->mSRSupport == 0){
                svs[i]->mPESupport = std::min(svs[i]->mPEWithPtn, svs[i]->mPESupport);
                sts->mSpnCnts[svs[i]->mID].mAltCnt = std::min(sts->mSpnCnts[svs[i]->mID].mAltCnt, svs[i]->mPESupport);
                sts->mTotalAltCnts[svs[i]->mID] =  std::min(sts->mSpnCnts[svs[i]->mID].mAltCnt, sts->mTotalAltCnts[svs[i]->mID]);
            }
        }
    }
    util::loginfo("Beg write final sv supporting bam");
    // write to result
    samFile* ifp = sam_open(mOpt->bamout.c_str(), "r");
    h = sam_hdr_read(ifp);
    b = bam_init1();
    while(sam_read1(ifp, h, b) >= 0){
        std::string qname = bamutil::getQName(b);
        auto ritr = rssm.find(qname);
        if(b->core.flag & BAM_FREAD1){
            bam_aux_update_int(b, "PH", ritr->second->mR1PHit);
            bam_aux_update_int(b, "SH", ritr->second->mR1SHit);
        }else{
            bam_aux_update_int(b, "PH", ritr->second->mR2PHit);
            bam_aux_update_int(b, "SH", ritr->second->mR2SHit);
        }
        auto iter = drec.find(qname);
        if(iter == drec.end()){
            assert(sam_write1(tmpSamFp, h, b) >= 0);
        }else{
            if(iter->second == 1 && (b->core.flag & BAM_FREAD2)){
                assert(sam_write1(tmpSamFp, h, b) >= 0);
            }else if(iter->second == 2 && (b->core.flag & BAM_FREAD1)){
                assert(sam_write1(tmpSamFp, h, b) >= 0);
            }
        }
    }
    sam_close(ifp);
    sam_close(tmpSamFp);
    // rename
    rename(tmpBam.c_str(), mOpt->bamout.c_str());
    // release mem
    bam_hdr_destroy(h);
    bam_destroy1(b);
    util::loginfo("End write final sv supporting bam");
}
