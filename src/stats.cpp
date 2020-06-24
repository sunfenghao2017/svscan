#include "stats.h"

Stats::Stats(Options* opt, int32_t n){
    mOpt = opt;
    init(n);
}

void Stats::init(int n){
    mJctCnts.resize(n, JunctionCount());
    mSpnCnts.resize(n, SpanningCount());
    mTotalAltCnts.resize(n, 0);
}

uint32_t Stats::getAlignmentQual(Matrix2D<char>* alnResult, const uint8_t* qual){
    int32_t baseQualSum = 0;
    int32_t readPos = 0;
    int32_t alignedBases = 0;
    for(int j = 0; j < alnResult->ncol(); ++j){
        if(alnResult->get(1, j) != '-'){
            if(alnResult->get(0, j) != '-'){
                ++alignedBases;
                baseQualSum += qual[readPos];
            }
            ++readPos;
        }
    }
    return baseQualSum/alignedBases;
}

bool Stats::validAlignment(Matrix2D<char>* alnResult, int32_t bppos, int32_t seqlen, int32_t minoff){
    int beg = 0, end = seqlen - 1;
    // find internal match range beg
    for(int i = 0; i < alnResult->ncol() - 3; ++i){
        if(alnResult->get(1, i) != '-') ++beg;
        if((alnResult->get(0, i) != '-') &&
           (alnResult->get(1, i) != '-') &&
           (alnResult->get(0, i+1) != '-') &&
           (alnResult->get(1, i+1) != '-') &&
           (alnResult->get(0, i+2) != '-') &&
           (alnResult->get(1, i+2) != '-')){
            break;
        }
    }
    // find internal match range end
    for(int i = alnResult->ncol() - 1; i >= 2; --i){
        if(alnResult->get(1, i) != '-') --end;
        if((alnResult->get(0, i) != '-') &&
           (alnResult->get(1, i) != '-') &&
           (alnResult->get(0, i-1) != '-') &&
           (alnResult->get(1, i-1) != '-') &&
           (alnResult->get(0, i-2) != '-') &&
           (alnResult->get(1, i-2) != '-')){
            break;
        }
    }
    // check validated bppos or not
    return bppos - beg > minoff && end - bppos > minoff;
}

Stats* Stats::merge(const std::vector<Stats*>& sts, int32_t n, Options* opt){
    Stats* ret = new Stats(opt, n);
    if(sts.size() == 0) return ret;
    for(int32_t j = 0; j < n; ++j){
        for(uint32_t i = 0; i < sts.size(); ++i){
            // Total
            ret->mTotalAltCnts[j] += sts[i]->mTotalAltCnts[j];
            // SC
            ret->mJctCnts[j].mFPIns += sts[i]->mJctCnts[j].mFPIns;
            ret->mJctCnts[j].mAltCnt += sts[i]->mJctCnts[j].mAltCnt;
            ret->mJctCnts[j].mRefCntBeg += sts[i]->mJctCnts[j].mRefCntBeg;
            ret->mJctCnts[j].mRefCntEnd += sts[i]->mJctCnts[j].mRefCntEnd;
            for(auto iter = sts[i]->mJctCnts[j].mAltQual.begin(); iter != sts[i]->mJctCnts[j].mAltQual.end(); ++iter){
                auto jter = ret->mJctCnts[j].mAltQual.find(iter->first);
                if(jter == ret->mJctCnts[j].mAltQual.end()){
                    ret->mJctCnts[j].mAltQual[iter->first] = iter->second;
                }else{
                    jter->second += iter->second;
                }
            }
            for(auto iter = sts[i]->mJctCnts[j].mRefQualBeg.begin(); iter != sts[i]->mJctCnts[j].mRefQualBeg.end(); ++iter){
                auto jter = ret->mJctCnts[j].mRefQualBeg.find(iter->first);
                if(jter == ret->mJctCnts[j].mRefQualBeg.end()){
                    ret->mJctCnts[j].mRefQualBeg[iter->first] = iter->second;
                }else{
                    jter->second += iter->second;
                }
            }
            for(auto iter = sts[i]->mJctCnts[j].mRefQualEnd.begin(); iter != sts[i]->mJctCnts[j].mRefQualEnd.end(); ++iter){
                auto jter = ret->mJctCnts[j].mRefQualEnd.find(iter->first);
                if(jter == ret->mJctCnts[j].mRefQualEnd.end()){
                    ret->mJctCnts[j].mRefQualEnd[iter->first] = iter->second;
                }else{
                    jter->second += iter->second;
                }
            }
            // DP
            ret->mSpnCnts[j].mAltCnt += sts[i]->mSpnCnts[j].mAltCnt;
            ret->mSpnCnts[j].mRefCntBeg += sts[i]->mSpnCnts[j].mRefCntBeg;
            ret->mSpnCnts[j].mRefCntEnd += sts[i]->mSpnCnts[j].mRefCntEnd;
            for(auto iter = sts[i]->mSpnCnts[j].mAltQual.begin(); iter != sts[i]->mSpnCnts[j].mAltQual.end(); ++iter){
                auto jter = ret->mSpnCnts[j].mAltQual.find(iter->first);
                if(jter == ret->mSpnCnts[j].mAltQual.end()){
                    ret->mSpnCnts[j].mAltQual[iter->first] = iter->second;
                }else{
                    jter->second += iter->second;
                }
            }
            for(auto iter = sts[i]->mSpnCnts[j].mRefQualBeg.begin(); iter != sts[i]->mSpnCnts[j].mRefQualBeg.end(); ++iter){
                auto jter = ret->mSpnCnts[j].mRefQualBeg.find(iter->first);
                if(jter == ret->mSpnCnts[j].mRefQualBeg.end()){
                    ret->mSpnCnts[j].mRefQualBeg[iter->first] = iter->second;
                }else{
                    jter->second += iter->second;
                }
            }
            for(auto iter = sts[i]->mSpnCnts[j].mRefQualEnd.begin(); iter != sts[i]->mSpnCnts[j].mRefQualEnd.end(); ++iter){
                auto jter = ret->mSpnCnts[j].mRefQualEnd.find(iter->first);
                if(jter == ret->mSpnCnts[j].mRefQualEnd.end()){
                    ret->mSpnCnts[j].mRefQualEnd[iter->first] = iter->second;
                }else{
                    jter->second += iter->second;
                }
            }
        }
    }
    return ret;
}

void Stats::stat(const SVSet& svs, const ContigBpRegions& bpRegs, const ContigSpanPoints& spPts, const RegItemCnt& regInfo, cgranges_t* ctgCgr){
    int32_t refIdx = regInfo.mTid;
    int32_t chrBeg = regInfo.mBeg;
    int32_t chrEnd = regInfo.mEnd;
    samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    util::loginfo("Beg gathering coverage information on contig: " + std::string(h->target_name[refIdx]) +
                  " [" + std::to_string(chrBeg) + "," + std::to_string(chrEnd) + "]", mOpt->logMtx);
    // Flag breakpoint regions
    std::set<bool> bpOccupied;
    for(uint32_t i = 0; i < bpRegs[refIdx].size(); ++i){
        for(int32_t k = bpRegs[refIdx][i].mRegStart; k < bpRegs[refIdx][i].mRegEnd; ++k) bpOccupied.insert(k);
    }
    // Flag spanning breakpoints
    std::set<bool> spanBp;
    for(uint32_t i = 0; i < spPts[refIdx].size(); ++i) spanBp.insert(spPts[refIdx][i].mBpPos);
    // Count reads
    hts_idx_t* idx = sam_index_load(fp, mOpt->bamfile.c_str());
    hts_itr_t* itr = sam_itr_queryi(idx, refIdx, chrBeg, chrEnd);
    bam1_t* b = bam_init1();
    AlignConfig alnCfg(5, -4, -4, -4, false, true);   
    uint16_t COV_STAT_SKIP_MASK = (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY);
    if(mOpt->rnamode) COV_STAT_SKIP_MASK |= BAM_FMUNMAP;
    auto itbp = bpRegs[refIdx].begin();
    auto itspnr = spPts[refIdx].begin();
    auto itspna = spPts[refIdx].begin();
    auto ltbp = itbp;
    auto ltsp = itspnr;
    auto ltap = itspna;
    int32_t lpos = -1;
    int32_t lsprp = -1;
    int32_t lspap = -1;
    while(sam_itr_next(fp, itr, b) >= 0){
#ifdef DEBUG
        bool hitdbgr = false;
        if(mOpt->debug & DEBUG_FRESR){
            hitdbgr = (bam_get_qname(b) == mOpt->qndbg);
            if(hitdbgr) std::cout << "target read visited" << std::endl;
        }
#endif
        if(regInfo.mInterleved && b->core.pos < chrBeg){
#ifdef DEBUG
            if(hitdbgr) std::cout << "target interleaved and skip" << std::endl;
#endif
            continue;
        }
        if(b->core.flag & COV_STAT_SKIP_MASK){
#ifdef DEBUG
            if(hitdbgr) std::cout << "target skiped by coverage stat skip MASK: " << b->core.flag << ":" << COV_STAT_SKIP_MASK << std::endl;
#endif
            continue;
        }
        if(b->core.qual < mOpt->filterOpt->mMinGenoQual){
#ifdef DEBUG
            if(hitdbgr) std::cout << "target skiped by mapQ check: " << b->core.qual << ":" << mOpt->filterOpt->mMinGenoQual << std::endl;
#endif
            continue;
        }
        if(!cr_isoverlap(ctgCgr, 
                         h->target_name[b->core.tid], 
                         std::max((BIGD_TYPE)0, b->core.pos - mOpt->libInfo->mMaxNormalISize), 
                         std::min((int32_t)(b->core.pos + mOpt->libInfo->mMaxNormalISize), (int32_t)h->target_len[b->core.tid]))){
#ifdef DEBUG
            if(hitdbgr){
                std::cout << "target read region overlap failed" << std::endl;
            }
#endif
           continue;
        }
#ifdef DEBUG
        if((mOpt->debug & DEBUG_FRESR) && hitdbgr) std::cout << "target read passed qc and region overlap check" << std::endl;
#endif
        // Count aligned basepair (small InDels)
        int32_t leadingSC = 0, leadingHC = 0;
        int32_t tailingSC = 0, tailingHC = 0;
        int32_t rp = b->core.pos; // reference pos
        int32_t ep = rp;
        uint32_t* cigar = bam_get_cigar(b);
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            int opint = bam_cigar_op(cigar[i]);
            int oplen = bam_cigar_oplen(cigar[i]);
            switch(opint){
                case BAM_CMATCH: case BAM_CDIFF: case BAM_CEQUAL: case BAM_CDEL: case BAM_CREF_SKIP:
                    ep += oplen;
                    break;
                case BAM_CSOFT_CLIP:
                    if(i == 0) leadingSC = oplen;
                    else tailingSC = oplen;
                    break;
                case BAM_CHARD_CLIP:
                    if(i == 0) leadingHC = oplen;
                    else tailingHC = oplen;
                    break;
                default:
                    break;
            }
        }
        if(leadingHC || tailingHC) continue; // skip reads with any hardclipings
        if(leadingSC && tailingSC) continue; // skip reads with both leading and tailing softclips
        int32_t bpapos = rp;
        if(tailingSC) bpapos = ep;
        bool assigned = false;
        uint8_t* sa = bam_aux_get(b, "SA");
        int32_t stid = -1, irpos = -1, erpos = -1, sscl = 0, sscr = 0, bpbpos = -1, scsvt = -1, seqmatch = 0;
        bool safwd = false, orifwd = !(b->core.flag & BAM_FREVERSE), sahdc = false;
        std::string sastr;
        if(sa){ // skip reads with cliped part in repeat regions
            sastr = bam_aux2Z(sa);
            std::vector<std::string> cvs;
            std::vector<std::string> vstr;
            util::split(sastr, cvs, ";");
            if(cvs[1].empty()){
                if(cvs[0].find_first_of("SH") == cvs[0].find_last_of("SH")) sastr = cvs[0];
                else sastr = "";
            }else{
                std::vector<int32_t> mvidx;
                for(uint32_t cvidx = 0; cvidx < cvs.size() - 1; ++cvidx){
                    util::split(cvs[cvidx], vstr, ",");
                    if(vstr[3].find_first_of("SH") == vstr[3].find_last_of("SH")){
                        mvidx.push_back(cvidx);
                    }
                }
                if(mvidx.size() == 1){
                    sastr = cvs[mvidx[0]];
                }else{
                    sastr = "";
                }
            }
            if(sastr.empty()){
                if(leadingSC + tailingSC < mOpt->filterOpt->mMinGoodSRLen){
                    bam_aux_del(b, sa);
                    sa = NULL;
                }else{
                    continue;
                }
            }else{
                util::split(sastr, vstr, ",");
                stid = bam_name2id(h, vstr[0].c_str());
                irpos = std::atoi(vstr[1].c_str()) - 1;
                erpos = irpos;
                safwd = (vstr[2][0] == '+');
                char* scg = const_cast<char*>(vstr[3].c_str());
                int32_t stotlen = 0;
                seqmatch = 0;
                while(*scg && *scg != '*'){
                    long num = 0;
                    if(std::isdigit((int)*scg)){
                        num = std::strtol(scg, &scg, 10);
                        stotlen += num;
                    }
                    switch(*scg){
                        case 'S':
                            if(stotlen == num) sscl = num;
                            else sscr = num;
                            break;
                        case 'H':
                            sahdc = true;
                            break;
                        case 'M': case '=': case 'X': case 'D': case 'N':
                            erpos += num;
                            seqmatch += num;
                            break;
                        default:
                            break;
                    }
                    ++scg;
                }
                bpbpos = irpos;
                if(sscr) bpbpos = erpos;
                if(seqmatch < mOpt->filterOpt->mMinGoodSRLen && seqmatch > (leadingSC + tailingSC)){
                    bam_aux_del(b, sa);
                    sa = NULL;
                }
                if(sa && stid == b->core.tid){
                    if(bpapos < bpbpos){
                        scsvt = svutil::getSRSASVT(bpapos, leadingSC, orifwd, bpbpos, sscl, safwd, true,
                                                   mOpt->filterOpt->mMaxReadSep, mOpt->filterOpt->mMinRefSep);
                    }else{
                        scsvt = svutil::getSRSASVT(bpbpos, sscl, safwd, bpapos, leadingSC, orifwd, true,
                                                   mOpt->filterOpt->mMaxReadSep, mOpt->filterOpt->mMinRefSep);
                    }
                }else{
                    if(b->core.tid < stid){
                        scsvt = svutil::getSRSASVT(bpapos, leadingSC, orifwd, bpbpos, sscl, safwd, false,
                                                   mOpt->filterOpt->mMaxReadSep, mOpt->filterOpt->mMinRefSep);
                    }else{
                        scsvt = svutil::getSRSASVT(bpbpos, sscl, safwd, bpapos, leadingSC, orifwd, false,
                                                   mOpt->filterOpt->mMaxReadSep, mOpt->filterOpt->mMinRefSep);
                    }
                }
            }
        }
        std::set<int32_t> sptids;
        int32_t svt = DPBamRecord::getSVType(b, mOpt);
        int32_t sst = DPBamRecord::getSVType(b);
        // Check read length for junction annotation
        if(b->core.l_qseq > 2 * mOpt->filterOpt->mMinFlankSize){
            bool bpvalid = false;
            int32_t rbegin = std::max((BIGD_TYPE)0, b->core.pos - leadingSC);
            int32_t rend = std::min(rbegin + b->core.l_qseq, (int32_t)h->target_len[refIdx]);
            for(int32_t k = rbegin; k < rend; ++k){
                if(bpOccupied.find(k) != bpOccupied.end()){
                    bpvalid = true;
                    break;
                }
            }
            if(bpvalid){
                bool onlySupportIns = true;
                std::vector<int32_t> supportInsID;
                std::map<int32_t, std::pair<float, bool>> supportSrsID;
                // get sequence
                std::string readOri;
                std::string readSeq;
                // Fetch all relevant SVs
                if(rbegin != lpos){
                    itbp = std::lower_bound(bpRegs[refIdx].begin(), bpRegs[refIdx].end(), BpRegion(rbegin));
                    ltbp = itbp;
                    lpos = rbegin;
                }else{// if same pos, do not search again, searching is really a time consuming work...
                    itbp = ltbp;
                }
                for(; itbp != bpRegs[refIdx].end() && rend >= itbp->mBpPos; ++itbp){
                    // Read spans breakpoint, if this read mapping range contains itbp->mBpPos ± mMinFlankSize
                    if(rbegin + mOpt->filterOpt->mMinFlankSize <= itbp->mBpPos && rend >= itbp->mBpPos + mOpt->filterOpt->mMinFlankSize){
                        if(!(leadingSC + tailingSC) && b->core.tid == b->core.mtid && sst >= 2){// REF Type, no realignment needed
                            bool add2ref = true;
                            if((b->core.mpos + mOpt->filterOpt->mMinFlankSize <= itbp->mBpPos) &&
                               (b->core.mpos + b->core.l_qseq >= itbp->mBpPos + mOpt->filterOpt->mMinFlankSize)){
                                if(b->core.flag & BAM_FREAD2) add2ref = true;
                                else add2ref = false;
                            }
                            if(add2ref && b->core.qual >= mOpt->filterOpt->mMinGenoQual){
                                mOpt->logMtx.lock();
                                if(itbp->mIsSVEnd){
                                    mJctCnts[itbp->mID].mRefCntEnd  += 1;
                                    if(mOpt->writebcf){
                                        auto qiter = mJctCnts[itbp->mID].mRefQualEnd.find(b->core.qual);
                                        if(qiter == mJctCnts[itbp->mID].mRefQualEnd.end()){
                                            mJctCnts[itbp->mID].mRefQualEnd[b->core.qual] = 1;
                                        }else{
                                            qiter->second += 1;
                                        }
                                    }
                                }else{
                                    mJctCnts[itbp->mID].mRefCntBeg += 1;
                                    if(mOpt->writebcf){
                                        auto qiter = mJctCnts[itbp->mID].mRefQualEnd.find(b->core.qual);
                                        if(qiter == mJctCnts[itbp->mID].mRefQualEnd.end()){
                                            mJctCnts[itbp->mID].mRefQualEnd[b->core.qual] = 1;
                                        }else{
                                            qiter->second += 1;
                                        }
                                    }
                                }
                                mOpt->logMtx.unlock();
                            }
                            continue;
                        }
                        if(leadingSC + tailingSC < mOpt->filterOpt->mMinRealnFlkLen) continue; // skip reads with too short softclips
                        bool seedgot = false;
                        int seedoff = 0;
                        if(sa){
#ifdef DEBUG
                            if((mOpt->debug & DEBUG_FRESR) &&  hitdbgr){
                                std::cout << "scsvt: " << scsvt << "\n";
                                std::cout << "itbp->mSVT: " << scsvt << "\n";
                                std::cout << "sahdc: " << sahdc << "\n";
                                std::cout << "itbp->mIsSVEnd: " << itbp->mIsSVEnd << "\n";
                                std::cout << "bpbpos: " << bpbpos << "\n";
                                std::cout << "svs[itbp->mID]->mSVStart: " << svs[itbp->mID]->mSVStart << "\n";
                                std::cout << "svs[itbp->mID]->mSVEnd: " << svs[itbp->mID]->mSVEnd << "\n";
                            }
#endif
                            if((scsvt != itbp->mSVT) || sahdc) continue;
                            if(itbp->mIsSVEnd){
                                if(std::abs(bpbpos - svs[itbp->mID]->mSVStart) >  mOpt->filterOpt->mMaxReadSep + svs[itbp->mID]->mHomLen || svs[itbp->mID]->mChr1 != stid){
                                    continue;
                                }
                                seedoff = std::abs(svs[itbp->mID]->mSVStart - bpbpos);
                            }else{
                                if(std::abs(bpbpos - svs[itbp->mID]->mSVEnd) > mOpt->filterOpt->mMaxReadSep + svs[itbp->mID]->mHomLen || svs[itbp->mID]->mChr2 != stid){
                                    continue;
                                }
                                seedoff = std::abs(svs[itbp->mID]->mSVEnd - bpbpos);
                            }
                            seedgot = true;
                        }
#ifdef DEBUG
                            if((mOpt->debug & DEBUG_FRESR) &&  hitdbgr){
                                std::cout << "seedgot: " << seedgot << "\n";
                                std::cout << "itbp->mSVT: " << scsvt << "\n";
                                std::cout << "seedoff: " << seedoff << "\n";
                            }
#endif
                        if(seedgot){
                            assigned = true;
                            if(itbp->mSVT == 4) supportInsID.push_back(itbp->mID);
                            if(itbp->mSVT != 4) onlySupportIns = false;
                            if(b->core.qual >= mOpt->filterOpt->mMinGenoQual) supportSrsID[itbp->mID] = {seedoff, itbp->mIsSVEnd};
#ifdef DEBUG
                            if((mOpt->debug & DEBUG_FRESR) &&  hitdbgr){
                                std::cout << "seed got good" << "\n";
                            }
#endif
                            continue;
                        }
                        if(svt >= 0 && svt != itbp->mSVT){
#ifdef DEBUG
                            if((mOpt->debug & DEBUG_FRESR) &&  hitdbgr){
                                std::cout << "non seed read with invalid svtype skipped" << "\n";
                                std::cout << "svt: " << svt << "\n";
                                std::cout << "bpsvt: " <<  svs[itbp->mID]->mSVT << "\n";
                            }
#endif

                            continue; // non-compatible svtype
                        }
                        // breakpoint check for rescue reads
                        if(itbp->mIsSVEnd){
                            if(std::abs(svs[itbp->mID]->mSVEnd - bpapos) >= mOpt->filterOpt->mMaxReadSep){
#ifdef DEBUG
                                if((mOpt->debug & DEBUG_FRESR) &&  hitdbgr){
                                    std::cout << "non seed read at sv end, too farwary skipped" << "\n";
                                    std::cout << "read svpos: " << bpapos << "\n";
                                    std::cout << "distance: " << svs[itbp->mID]->mSVEnd - bpapos << "\n";
                                }
#endif
                                continue;
                            }
                        }else{
                            if(std::abs(svs[itbp->mID]->mSVStart - bpapos) >= mOpt->filterOpt->mMaxReadSep){
#ifdef DEBUG
                                if((mOpt->debug & DEBUG_FRESR) &&  hitdbgr){
                                    std::cout << "non seed read at sv beg, too farwary skipped" << "\n";
                                    std::cout << "read svpos: " << bpapos << "\n";
                                    std::cout << "distance: " << svs[itbp->mID]->mSVStart << "\n";
                                }
#endif
                                continue;
                            }
                        }
                        // possible ALT
                        std::string consProbe = itbp->mIsSVEnd ? svs[itbp->mID]->mProbeEndC : svs[itbp->mID]->mProbeBegC;
                        if(readOri.empty()) readOri = bamutil::getSeq(b); // then fetch read to do realign
                        readSeq = readOri;
                        std::string adjscseq;
                        bool recov = SRBamRecord::adjustOrientation(readSeq, itbp->mIsSVEnd, itbp->mSVT);
                        int adjbppos = 0;
                        if(leadingSC){
                            if(recov){
                                adjbppos = readSeq.length() - leadingSC;
                                adjscseq = readSeq.substr(adjbppos);
                            }else{
                                adjbppos = leadingSC;
                                adjscseq = readSeq.substr(0, leadingSC);
                            }
                        }else{
                            if(recov){
                                adjbppos = tailingSC;
                                adjscseq = readSeq.substr(0, tailingSC);
                            }else{
                                adjbppos = readSeq.length() - tailingSC;
                                adjscseq = readSeq.substr(adjbppos);
                            }
                        }
                        // Compute alignment to alternative haplotype
                        Aligner* altAligner = new Aligner(consProbe, readSeq, &alnCfg);
                        Matrix2D<char>* altResult = new Matrix2D<char>();
                        int alnScore = altAligner->needle(altResult);
                        int matchThreshold = mOpt->filterOpt->mFlankQuality * consProbe.size() * alnCfg.mMatch + (1 - mOpt->filterOpt->mFlankQuality) * consProbe.size() * alnCfg.mMisMatch;
                        double scoreAlt = (double)alnScore / (double)matchThreshold;
                        // check match range on read
                        if(scoreAlt >= mOpt->filterOpt->mMinSRResScore){
                            if(!validAlignment(altResult, adjbppos, readSeq.length())){
#ifdef DEBUG
                                if((mOpt->debug & DEBUG_FRESR) &&  hitdbgr){
                                    std::cout << "non seed read invalid realignment on cns" << "\n";
                                    std::cout << *altResult << "\n";
                                    std::cout << "adjbppos: " << adjbppos << "\n";
                                    std::cout << "distance: " << svs[itbp->mID]->mSVStart << "\n";
                                }
#endif
                                // free resouces and go to next round
                                delete altAligner; altAligner = NULL; delete altResult; altResult = NULL;
                                continue;
                            }else{
#ifdef DEBUG
                                if((mOpt->debug & DEBUG_FRESR) &&  hitdbgr){
                                    std::cout << "non seed read alignment good" << "\n";
                                    std::cout << "scoreAlt: " << scoreAlt << "\n";
                                    std::cout << "mathThr: " << matchThreshold << "\n";
                                    std::cout << "scoreThr: " << mOpt->filterOpt->mMinSRResScore << "\n";
                                    std::cout << "svsrpt: " << "\n";
                                    std::cout << svs[itbp->mID] << std::endl;
                                }
#endif
                                // just free resources
                                delete altAligner; altAligner = NULL; delete altResult; altResult = NULL;
                                Aligner* ascAligner = new Aligner(adjscseq, svs[itbp->mID]->mConsensus, &alnCfg);
                                Matrix2D<char>* ascResult = new Matrix2D<char>();
                                double ascScore = ascAligner->needle(ascResult);
                                double ascMThre = mOpt->filterOpt->mFlankQuality * adjscseq.size() * alnCfg.mMatch + (1 - mOpt->filterOpt->mFlankQuality) * adjscseq.size() * alnCfg.mMisMatch;
                                double ascSAlt = ascScore/ascMThre;
                                if(ascSAlt < mOpt->filterOpt->mMinSRResScore){
#ifdef DEBUG
                                    if((mOpt->debug & DEBUG_FRESR) &&  hitdbgr){
                                        std::cout << "non seed read alignment to flank sequence bad" << "\n";
                                        std::cout << "ascScore: " << ascScore << "\n";
                                        std::cout << "ascMThre: " << ascMThre << "\n";
                                        std::cout << "scoreThr: " << mOpt->filterOpt->mMinSRResScore << "\n";
                                        std::cout << *ascResult << std::endl;
                                        std::cout << svs[itbp->mID] << std::endl;
                                    }
#endif
                                    scoreAlt = 0.0;
                                }
                                delete ascAligner; ascAligner = NULL; delete ascResult; ascResult = NULL;
                            }
                        }else{
                            // free previous resources
                            delete altAligner; altAligner = NULL; delete altResult; altResult = NULL;
                            // Any confident alignment?
                            if(scoreAlt < mOpt->filterOpt->mMinSRResScore &&
                               svs[itbp->mID]->mProbeEndA.size() && 
                               (leadingSC + tailingSC > (int32_t)(svs[itbp->mID]->mBpInsSeq.size() + mOpt->filterOpt->mMinInsFlkLen))){
                                consProbe = itbp->mIsSVEnd ? svs[itbp->mID]->mProbeEndA : svs[itbp->mID]->mProbeBegA;
                                Aligner* secAligner = new Aligner(consProbe, readSeq, &alnCfg);
                                Matrix2D<char>* secResult = new Matrix2D<char>();
                                alnScore = secAligner->needle(secResult);
                                matchThreshold = mOpt->filterOpt->mFlankQuality * consProbe.size() * alnCfg.mMatch + (1 - mOpt->filterOpt->mFlankQuality) * consProbe.size() * alnCfg.mMisMatch;
                                scoreAlt = (double)alnScore / (double)matchThreshold;
                                // check match range on read
                                if(scoreAlt >= mOpt->filterOpt->mMinSRResScore){
                                    if(!validAlignment(secResult, adjbppos, readSeq.length())){
                                        // free resouces and go to next round 
                                        delete secAligner; secResult = NULL; delete secResult; secResult = NULL;
                                        continue;
                                    }else{
                                        // just free resources
                                        delete secAligner; secResult = NULL; delete secResult; secResult = NULL;
                                        Aligner* ascAligner = new Aligner(adjscseq, svs[itbp->mID]->mConsensus, &alnCfg);
                                        Matrix2D<char>* ascResult = new Matrix2D<char>();
                                        double ascScore = ascAligner->needle(ascResult);
                                        double ascMThre = mOpt->filterOpt->mFlankQuality * adjscseq.size() * alnCfg.mMatch + (1 - mOpt->filterOpt->mFlankQuality) * adjscseq.size() * alnCfg.mMisMatch;
                                        double ascSAlt = ascScore/ascMThre;
                                        if(ascSAlt < mOpt->filterOpt->mMinSRResScore) scoreAlt = 0.0;
                                        delete ascAligner; ascAligner = NULL; delete ascResult; ascResult = NULL;
                                    }
                                }else{
                                    delete secAligner; secResult = NULL; delete secResult; secResult = NULL;
                                }
                            }
                        }
                        if(scoreAlt >= mOpt->filterOpt->mMinSRResScore){
                            assigned = true;
                            if(itbp->mSVT == 4) supportInsID.push_back(itbp->mID);
                            if(itbp->mSVT != 4) onlySupportIns = false;
                            if(b->core.qual >= mOpt->filterOpt->mMinGenoQual) supportSrsID[itbp->mID] = {1.0/scoreAlt, itbp->mIsSVEnd};
                        }
                    }
                }
                if(supportSrsID.size()){
                    auto iit = supportSrsID.begin();
                    int32_t ixmx = iit->first;
                    float mmxhv = iit->second.first;
                    int32_t ixmcnt = 1;
                    if(supportSrsID.size() > 1){
                        // get the one with most confident hit
                        for(; iit != supportSrsID.end(); ++iit){
                            if(iit->second.first < mmxhv){
                                ixmx = iit->first;
                                mmxhv = iit->second.first;
                                ixmcnt = 1;
                            }else if(iit->second.first == mmxhv){
                                ixmcnt += 1;
                            }
                        }
                        // fix tie
                        if(ixmcnt > 1){
                            // din non-repeat region svs
                            std::vector<int32_t> nrps;
                            for(iit = supportSrsID.begin(); iit != supportSrsID.end(); ++iit){
                                if(iit->second.first == mmxhv &&
                                   svs[iit->first]->mRealnRet >= 0 && 
                                   svs[iit->first]->mRealnRet <= mOpt->fuseOpt->mWhiteFilter.mMaxRepHit){
                                    nrps.push_back(iit->first);
                                }
                            }
                            if(nrps.size() == 1) ixmx = nrps[0]; // only one non-repeat region hit
                            else{
                                if(nrps.size() == 0){ // all in repeat region
                                    // first run
                                    for(iit = supportSrsID.begin(); iit != supportSrsID.end(); ++iit){
                                        if(svs[iit->first]->mSRSupport > svs[ixmx]->mSRSupport){
                                            ixmx = iit->first;
                                        }
                                    }
                                    // second run
                                    for(iit = supportSrsID.begin(); iit != supportSrsID.end(); ++iit){
                                        if(svs[iit->first]->mSRSupport == svs[ixmx]->mSRSupport && 
                                           svs[iit->first]->mPESupport > svs[ixmx]->mPESupport){
                                           ixmx = iit->first;
                                        }
                                    }
                                }else{ // more than 2 non-repeat region
                                    ixmx = nrps[0];
                                    // first run
                                    for(uint32_t xxvid = 1; xxvid < nrps.size(); ++xxvid){
                                        if(svs[nrps[xxvid]]->mSRSupport > svs[ixmx]->mSRSupport){
                                            ixmx = nrps[xxvid];
                                        }
                                    }
                                    // second run
                                    for(uint32_t xxvid = 0; xxvid < nrps.size(); ++xxvid){
                                        if(svs[nrps[xxvid]]->mSRSupport == svs[ixmx]->mSRSupport &&
                                           svs[nrps[xxvid]]->mPESupport > svs[ixmx]->mPESupport){
                                            ixmx = nrps[xxvid];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // output
                    mOpt->logMtx.lock();
                    mJctCnts[ixmx].mAltCnt += 1;
                    mTotalAltCnts[ixmx] += 1;
                    if(mOpt->writebcf){
                        auto qiter = mJctCnts[ixmx].mAltQual.find(b->core.qual);
                        if(qiter == mJctCnts[ixmx].mAltQual.end()){
                            mJctCnts[ixmx].mAltQual[b->core.qual] = 1;
                        }else{
                            qiter->second += 1;
                        }
                    }
                    if(mOpt->fbamout){
                        bam_aux_update_int(b, "ZF", ixmx);
                        bam_aux_update_int(b, "ST", 0);
                        assert(sam_write1(mOpt->fbamout, h, b) >= 0);
                        sptids.insert(ixmx);
                    }
                    mOpt->logMtx.unlock();
                }
                if(supportInsID.size() > 0 && (!onlySupportIns)){
                    for(auto& insid : supportInsID) ++mJctCnts[insid].mFPIns;
                }
            }
        }
        // Read-count and spanning annotation
        if(sa || (!(b->core.flag & BAM_FPAIRED))) continue;
        if(leadingSC + tailingSC >= mOpt->filterOpt->mMinGoodSRLen) continue; // rubbish sequence in tail
        if(b->core.flag & BAM_FMUNMAP) continue; // mate unmapped, skip
        if(b->core.tid > b->core.mtid || (b->core.tid == b->core.mtid && b->core.pos > b->core.mpos)){// Second read in pair
            if(b->core.qual < mOpt->filterOpt->mMinGenoQual) continue; // Low quality pair
            // Spanning counting
            int32_t outerISize = b->core.pos + b->core.l_qseq - b->core.mpos;
            // Normal spanning pair
            if(sst == 2 && outerISize >= mOpt->libInfo->mMinNormalISize &&
               outerISize <= mOpt->libInfo->mMaxNormalISize && b->core.mtid == b->core.tid){
                // Take 80% of the outersize as the spanned interval
                int32_t spanlen = 0.8 * outerISize;
                int32_t pbegin = b->core.mpos;
                int32_t st = pbegin + 0.1 * outerISize;
                bool spanvalid = false;
                for(int32_t i = st; i < st + spanlen && i < (int32_t)h->target_len[refIdx]; ++i){
                    if(spanBp.find(i) != spanBp.end()){
                        spanvalid = true;
                        break;
                    }
                }
                if(spanvalid){
                    // Fetch all relevant SVs
                    if(st != lsprp){
                        itspnr = std::lower_bound(spPts[refIdx].begin(), spPts[refIdx].end(), SpanPoint(st));
                        lsprp = st;
                        ltsp = itspnr;
                    }else{ // not search again if st is the same...
                        itspnr = ltsp;
                    }
                    for(; itspnr != spPts[refIdx].end() && (st + spanlen) >= itspnr->mBpPos; ++itspnr){
                        mOpt->logMtx.lock();
                        if(itspnr->mIsSVEnd){
                            ++mSpnCnts[itspnr->mID].mRefCntEnd;
                            if(mOpt->writebcf){
                                auto qiter = mSpnCnts[itspnr->mID].mRefQualEnd.find(b->core.qual);
                                if(qiter == mSpnCnts[itspnr->mID].mRefQualEnd.end()){
                                    mSpnCnts[itspnr->mID].mRefQualEnd[b->core.qual] = 1;
                                }else{
                                    qiter->second += 1;
                                }
                            }
                        }else{
                            ++mSpnCnts[itspnr->mID].mRefCntBeg;
                            if(mOpt->writebcf){
                                auto qiter = mSpnCnts[itspnr->mID].mRefQualBeg.find(b->core.qual);
                                if(qiter == mSpnCnts[itspnr->mID].mRefQualBeg.end()){
                                    mSpnCnts[itspnr->mID].mRefQualBeg[b->core.qual] = 1;
                                }else{
                                    qiter->second += 1;
                                }
                            }
                        }
                        mOpt->logMtx.unlock();
                    }
                }
            }
            // Abnormal spanning coverage
            if(svt == -1) continue;
            // Spanning a breakpoint?
            bool spanvalid = false;
            DPBamRecord dpr(b, svt);
            int32_t svstart = -1, svend = -1, pbegin = -1, pend = -1;
            dpr.initClique(svstart, svend, svt);
            if(dpr.mCurTid == dpr.mMateTid){
                pbegin = svstart - mOpt->libInfo->mVarisize;
                pend = svend + mOpt->libInfo->mVarisize;
            }else{
                pbegin = svstart - mOpt->libInfo->mVarisize;
                pend = svstart + mOpt->libInfo->mVarisize;
            }
            for(int32_t i = pbegin; i < pend; ++i){
                if(spanBp.find(i) != spanBp.end()){
                    spanvalid = true;
                    break;
                }
            }
            if(spanvalid){
                std::map<int32_t, std::pair<int32_t, bool>> supportSpnID;
                if(lspap != pbegin){
                    itspna = std::lower_bound(spPts[refIdx].begin(), spPts[refIdx].end(), SpanPoint(pbegin));
                    lspap = pbegin;
                    ltap = itspna;
                }else{
                    itspna = ltap;
                }
                for(; itspna != spPts[refIdx].end() && pend >= itspna->mBpPos; ++itspna){
                    if(svt == itspna->mSVT && svs[itspna->mID]->mChr1 == b->core.tid && svs[itspna->mID]->mChr2 == b->core.mtid){
                        // valid dp in pe-mode
                        bool validDPE = true;
                        if(std::abs(svend - svs[itspna->mID]->mSVEnd) > mOpt->libInfo->mVarisize){
                            validDPE = false;
                        }
                        if(std::abs(svstart - svs[itspna->mID]->mSVStart) > mOpt->libInfo->mVarisize){
                            validDPE = false;
                        }
                        if(svs[itspna->mID]->mConsensus.size() && (leadingSC + tailingSC)){
                            if(itspna->mIsSVEnd){
                                if(std::abs(bpapos - svs[itspna->mID]->mSVEnd) > mOpt->filterOpt->mMaxReadSep){
                                    validDPE = false;
                                }
                            }else{
                               if(std::abs(bpapos - svs[itspna->mID]->mSVStart) > mOpt->filterOpt->mMaxReadSep){
                                   validDPE = false;
                               }
                            }
                        }
                        if(validDPE){
                            supportSpnID[itspna->mID] = {std::abs(svend - svs[itspna->mID]->mSVEnd) + std::abs(svstart - svs[itspna->mID]->mSVStart),itspna->mIsSVEnd};
                        }
                    }
                }
                if(supportSpnID.size()){
                    auto iit = supportSpnID.begin();
                    int32_t ixmx = iit->first;
                    int32_t mmxhv = iit->second.first;
                    int32_t ixmcnt = 1;
                    if(supportSpnID.size() > 1){
                        // get the one with most cnofident hit
                        for(; iit != supportSpnID.end(); ++iit){
                            if(iit->second.first < mmxhv){
                                ixmx = iit->first;
                                mmxhv = iit->second.first;
                                ixmcnt = 1;
                            }else if(iit->second.first == mmxhv){
                                ixmcnt += 1;
                            }
                        }
                        // fix tie
                        if(ixmcnt > 1){
                            // in non-repeat region svs
                            std::vector<int32_t> nrps;
                            for(iit = supportSpnID.begin(); iit != supportSpnID.end(); ++iit){
                                if(iit->second.first == mmxhv &&
                                   svs[iit->first]->mRealnRet >= 0 && 
                                   svs[iit->first]->mRealnRet <= mOpt->fuseOpt->mWhiteFilter.mMaxRepHit){
                                    nrps.push_back(iit->first);
                                }
                            }
                            if(nrps.size() == 1) ixmx = nrps[0]; // only one non-repeat region hit
                            else{
                                if(nrps.size() == 0){ // all in repeat region
                                    // first run
                                    for(iit = supportSpnID.begin(); iit != supportSpnID.end(); ++iit){
                                        if(svs[iit->first]->mSRSupport > svs[ixmx]->mSRSupport){
                                            ixmx = iit->first;
                                        }
                                    }
                                    // second run
                                    for(iit = supportSpnID.begin(); iit != supportSpnID.end(); ++iit){
                                        if(svs[iit->first]->mSRSupport == svs[ixmx]->mSRSupport && 
                                           svs[iit->first]->mPESupport > svs[ixmx]->mPESupport){
                                           ixmx = iit->first;
                                        }
                                    }
                                }else{ // more than 2 non-repeat region
                                    ixmx = nrps[0];
                                    // first run
                                    for(uint32_t xxvid = 1; xxvid < nrps.size(); ++xxvid){
                                        if(svs[nrps[xxvid]]->mSRSupport > svs[ixmx]->mSRSupport){
                                            ixmx = nrps[xxvid];
                                        }
                                    }
                                    // second run
                                    for(uint32_t xxvid = 0; xxvid < nrps.size(); ++xxvid){
                                        if(svs[nrps[xxvid]]->mSRSupport == svs[ixmx]->mSRSupport &&
                                           svs[nrps[xxvid]]->mPESupport > svs[ixmx]->mPESupport){
                                            ixmx = nrps[xxvid];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // output
                    mOpt->logMtx.lock();
                    mSpnCnts[ixmx].mAltCnt += 1;
                    if(!assigned) mTotalAltCnts[ixmx] += 1;
                    if(mOpt->writebcf){
                        auto qiter = mSpnCnts[ixmx].mAltQual.find(b->core.qual);
                        if(qiter == mSpnCnts[ixmx].mAltQual.end()){
                            mSpnCnts[ixmx].mAltQual[b->core.qual] = 1;
                        }else{
                            qiter->second += 1;
                        }
                    }
                    if((!assigned) && mOpt->fbamout){
                        bam_aux_update_int(b, "ZF", ixmx);
                        bam_aux_update_int(b, "ST", 1);
                        assert(sam_write1(mOpt->fbamout, h, b) >= 0);
                        sptids.insert(ixmx);
                    }
                    mOpt->logMtx.unlock();
                }
            }
        }
    }
    util::loginfo("End gathering coverage information on contig: " + std::string(h->target_name[refIdx]) +
                  " [" + std::to_string(chrBeg) + "," + std::to_string(chrEnd) + "]", mOpt->logMtx);
    // Clean-up
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
    if(itr) hts_itr_destroy(itr);
    hts_idx_destroy(idx);
}
