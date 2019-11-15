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

Stats* Annotator::covAnnotate(std::vector<SVRecord>& svs){
    // Open file handler
    samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
    hts_set_fai_filename(fp, mOpt->genome.c_str());
    bam_hdr_t* h = sam_hdr_read(fp);
    faidx_t* fai = fai_load(mOpt->genome.c_str());
    hts_idx_t* idx = sam_index_load(fp, mOpt->bamfile.c_str());
    // Find Ns in reference genome
    util::loginfo("Beg gathering N regions in reference");
    std::vector<std::set<std::pair<int32_t, int32_t>>> nreg(h->n_targets);
    for(auto& refIndex : mOpt->svRefID){
        int32_t refSeqLen = -1;
        char* seq = faidx_fetch_seq(fai, h->target_name[refIndex], 0, h->target_len[refIndex], &refSeqLen);
        bool nrun = false;
        int nstart = refSeqLen;
        for(int32_t i = 0; i < refSeqLen; ++i){
            if(seq[i] != 'N' || seq[i] != 'n'){
                if(nrun){
                    nreg[refIndex].insert(std::make_pair(nstart, i - 1));
                    nrun = false;
                }
            }else{
                if(!nrun){
                    nstart = i;
                    nrun = true;
                }
            }
        }
        // Insert last possible Ns region
        if(nrun) nreg[refIndex].insert(std::make_pair(nstart, refSeqLen - 1));
        if(seq) free(seq);
    }
    util::loginfo("End gathering N regions in reference");
    util::loginfo("Beg extracting left/middle/right regions for each SV");
    // Add control regions
    std::vector<std::vector<CovRecord>> covRecs(h->n_targets); //coverage records of 3-part of each SV events
    int32_t lastID = svs.size();
    for(auto itsv = svs.begin(); itsv != svs.end(); ++itsv){
        int halfsize = (itsv->mSVEnd - itsv->mSVStart) / 2; 
        if(itsv->mSVT >= 4) halfsize = 500;
        // Left control region
        CovRecord covLeft;
        covLeft.mID = lastID + itsv->mID;
        covLeft.mStart = std::max(0, itsv->mSVStart - halfsize);
        covLeft.mEnd = itsv->mSVStart;
        auto itp = getFirstOverlap(nreg[itsv->mChr1], {covLeft.mStart, covLeft.mEnd});
        while(itp != nreg[itsv->mChr1].end()){
            covLeft.mStart = std::max(itp->first - halfsize, 0);
            covLeft.mEnd = itp->first;
            itp = getFirstOverlap(nreg[itsv->mChr1], {covLeft.mStart, covLeft.mEnd});
        }
        covRecs[itsv->mChr1].push_back(covLeft);
        // Actual SV region
        CovRecord covMiddle;
        covMiddle.mID = itsv->mID;
        covMiddle.mStart = itsv->mSVStart;
        covMiddle.mEnd = itsv->mSVEnd;
        if(itsv->mSVT >= 4){
            covMiddle.mStart = std::max(itsv->mSVStart - halfsize, 0);
            covMiddle.mEnd = std::min((int32_t)h->target_len[itsv->mChr2], itsv->mSVEnd + halfsize);
        }
        itsv->mSize = itsv->mSVEnd - itsv->mSVStart;
        if(itsv->mSVT == 4) itsv->mSize = itsv->mInsSeq.size();
        if(itsv->mSVT >= 5) itsv->mSize = 666666666;
        covRecs[itsv->mChr1].push_back(covMiddle);
        // Right control region
        CovRecord covRight;
        covRight.mID = 2 * lastID + itsv->mID;
        covRight.mStart = itsv->mSVEnd;
        covRight.mEnd = std::min((int32_t)h->target_len[itsv->mChr2], itsv->mSVEnd + halfsize);
        itp = getFirstOverlap(nreg[itsv->mChr2], {covRight.mStart, covRight.mEnd});
        while(itp != nreg[itsv->mChr2].end()){
            covRight.mStart = itp->second;
            covRight.mEnd = itp->second + halfsize;
            itp = getFirstOverlap(nreg[itsv->mChr2], {covRight.mStart, covRight.mEnd});
        }
        covRecs[itsv->mChr2].push_back(covRight);
    }
    // Sort Coverage Records
    for(auto& refIndex : mOpt->svRefID) std::sort(covRecs[refIndex].begin(), covRecs[refIndex].end());
    util::loginfo("End extracting left/middle/right regions for each SV");
    // Preprocess REF and ALT
    ContigBpRegions bpRegion(h->n_targets);
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
                regEnd = std::min(itsv->mSVEnd + mOpt->filterOpt->mMinFlankSize, (int32_t)h->target_len[itsv->mChr2]);
                bpPos = itsv->mSVEnd;
            }else{// SV starting position region
                regChr = itsv->mChr1;
                regStart = std::max(0, itsv->mSVStart - mOpt->filterOpt->mMinFlankSize);
                regEnd = std::min(itsv->mSVStart + mOpt->filterOpt->mMinFlankSize, (int32_t)h->target_len[itsv->mChr1]);
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
    // Clean-up
    sam_close(fp);
    hts_idx_destroy(idx);
    fai_destroy(fai);
    // Get coverage from each contig in parallel
    std::vector<Stats*> covStats(mOpt->svRefID.size(), NULL);
    std::vector<std::future<void>> statRets(mOpt->svRefID.size());
    int32_t i = 0;
    for(auto& refidx: mOpt->svRefID){
        covStats[i] = new Stats(mOpt, svs.size(), refidx);
        statRets[i] = mOpt->pool->enqueue(&Stats::stat, covStats[i], std::ref(svs), std::ref(covRecs), std::ref(bpRegion), std::ref(spanPoint));
        ++i;
    }
    for(auto& e: statRets) e.get();
    Stats* finalStat = Stats::merge(covStats, svs.size(), mOpt);
    for(auto& e: covStats) delete e;
    return finalStat;
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
    hts_itr_destroy(itr);
    if(rec.s) free(rec.s);
}

void Annotator::geneAnnoDNA(SVSet& svs, GeneInfoList& gl){
    gl.resize(svs.size());
    // split range list 
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    int32_t totalSV = svs.size();
    int32_t eachTSV = totalSV / mOpt->nthread;
    for(int32_t i = 0; i < mOpt->nthread; ++i){
        std::pair<int32_t, int32_t> p;
        p.first = i * eachTSV;
        p.second = (i + 1) * eachTSV;
        if(p.second <= totalSV) vpidx.push_back(p);
        else break;
    }
    if(vpidx.size()) vpidx[vpidx.size() - 1].second = svs.size();
    else vpidx.push_back({0, vpidx.size()});
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
    tbx_destroy(tbx);
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
    hts_itr_destroy(itr);
    if(rec.s) free(rec.s);
}

void Annotator::geneAnnoRNA(SVSet& svs, GeneInfoList& gl){
    gl.resize(svs.size());
    // split range list 
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    int32_t totalSV = svs.size();
    int32_t eachTSV = totalSV / mOpt->nthread;
    for(int32_t i = 0; i < mOpt->nthread; ++i){
        std::pair<int32_t, int32_t> p;
        p.first = i * eachTSV;
        p.second = (i + 1) * eachTSV;
        if(p.second <= totalSV) vpidx.push_back(p);
        else break;
    }
    if(vpidx.size()) vpidx[vpidx.size() - 1].second = svs.size();
    else vpidx.push_back({0, svs.size()});
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
    tbx_destroy(tbx);
    hts_close(fp);
}
