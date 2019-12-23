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
    const uint16_t COV_STAT_SKIP_MASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FMUNMAP);
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
        if(regInfo.mInterleved && b->core.pos < chrBeg) continue;
        if(b->core.flag & COV_STAT_SKIP_MASK) continue;
        if(b->core.qual < mOpt->filterOpt->mMinGenoQual) continue;
        if(!cr_isoverlap(ctgCgr, 
                         h->target_name[b->core.tid], 
                         std::max((BIGD_TYPE)0, b->core.pos - mOpt->libInfo->mMaxNormalISize), 
                         std::min((int32_t)(b->core.pos + mOpt->libInfo->mMaxNormalISize), (int32_t)h->target_len[b->core.tid]))){
           continue;
        }
        // Count aligned basepair (small InDels)
        int32_t leadingSC = 0, leadingHC = 0;
        int32_t tailingSC = 0, tailingHC = 0;
        int32_t rp = 0; // reference pos
        uint32_t* cigar = bam_get_cigar(b);
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            int opint = bam_cigar_op(cigar[i]);
            int oplen = bam_cigar_oplen(cigar[i]);
            if(opint == BAM_CMATCH || opint == BAM_CDIFF || opint == BAM_CEQUAL){
                rp += oplen;
            }else if(opint == BAM_CDEL){
                rp += oplen;
            }else if(opint == BAM_CREF_SKIP){
                rp += oplen;
            }else if(opint == BAM_CSOFT_CLIP){
                if(i == 0) leadingSC = oplen;
                else tailingSC = opint;
            }else if(opint == BAM_CHARD_CLIP){
                if(i == 0) leadingHC = oplen;
                else tailingHC = oplen;
            }
        }
        if(leadingHC || tailingHC) continue; // skip reads with any hardclipings
        if(leadingSC && tailingSC) continue; // skip reads with both leading and tailing softclips
        bool assigned = false;
        uint8_t* sa = bam_aux_get(b, "SA");
        std::string sastr;
        if(sa){ // skip reads with cliped part in repeat regions
            sastr = bam_aux2Z(sa);
            std::vector<std::string> cvs;
            std::vector<std::string> vstr;
            util::split(sastr, cvs, ";");
            if(cvs[1].empty()) sastr = cvs[0];
            else{
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
            if(sastr.empty()) continue;
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
                std::map<int32_t, bool> supportSrsID;
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
                        // possible ALT
                        std::string consProbe = itbp->mIsSVEnd ? svs[itbp->mID].mProbeEndC : svs[itbp->mID].mProbeBegC;
                        if(readOri.empty()) readOri = bamutil::getSeq(b); // then fetch read to do realign
                        readSeq = readOri;
                        SRBamRecord::adjustOrientation(readSeq, itbp->mIsSVEnd, itbp->mSVT);
                        // Compute alignment to alternative haplotype
                        Aligner* altAligner = new Aligner(consProbe, readSeq, &alnCfg);
                        Matrix2D<char>* altResult = new Matrix2D<char>();
                        int alnScore = altAligner->needle(altResult);
                        int matchThreshold = mOpt->filterOpt->mFlankQuality * consProbe.size() * alnCfg.mMatch + (1 - mOpt->filterOpt->mFlankQuality) * consProbe.size() * alnCfg.mMisMatch;
                        double scoreAlt = (double)alnScore / (double)matchThreshold;
                        // Any confident alignment?
                        if(scoreAlt > mOpt->filterOpt->mMinSRResScore){
                            bool validRSR = true;
                            if(sa){
                                std::vector<std::string> vstr;
                                util::split(sastr, vstr, ",");
                                int32_t stid = bam_name2id(h, vstr[0].c_str());
                                int32_t irpos = std::atoi(vstr[1].c_str());
                                int32_t erpos = irpos;
                                bool safwd = (vstr[2][0] == '+');
                                bool orifwd = !(b->core.flag & BAM_FREVERSE);
                                int32_t catt = itbp->mSVT;
                                if(itbp->mSVT >= 5) catt -= 5;
                                if(catt <= 1){
                                    if(safwd == orifwd) validRSR = false;
                                }
                                if(catt >= 2){
                                    if(safwd != orifwd) validRSR = false;
                                }
                                if(validRSR){
                                    char* scg = const_cast<char*>(vstr[3].c_str());
                                    int32_t sscl = 0, sscr = 0;
                                    int32_t stotlen = 0;
                                    while(validRSR && *scg && *scg != '*'){
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
                                                validRSR = false;
                                                break;
                                            case 'M': case '=': case 'X': case 'D': case 'N':
                                                erpos += num;
                                                break;
                                            default:
                                                break;
                                        }
                                        ++scg;
                                    }
                                    if(validRSR && ((sscl > 0) ^ (sscr > 0))){
                                        int32_t bppos = irpos;
                                        if(sscr) bppos = erpos;
                                        if(itbp->mIsSVEnd){
                                            if(bppos < svs[itbp->mID].mSVStart - mOpt->libInfo->mReadLen ||
                                               bppos > svs[itbp->mID].mSVStart + mOpt->libInfo->mReadLen ||
                                               svs[itbp->mID].mChr1 != stid){
                                                validRSR = false;
                                            }
                                        }else{
                                            if(bppos < svs[itbp->mID].mSVEnd - mOpt->libInfo->mReadLen ||
                                               bppos > svs[itbp->mID].mSVEnd + mOpt->libInfo->mReadLen ||
                                               svs[itbp->mID].mChr2 != stid){
                                                validRSR = false;
                                            }
                                        }
                                    }else{
                                        validRSR = false;
                                    }
                                }
                            }
                            if(validRSR){
                                assigned = true;
                                if(itbp->mSVT == 4) supportInsID.push_back(itbp->mID);
                                if(itbp->mSVT != 4) onlySupportIns = false;
                                if(b->core.qual >= mOpt->filterOpt->mMinGenoQual) supportSrsID[itbp->mID] = itbp->mIsSVEnd;
                            }
                        }
                        delete altAligner;
                        altAligner = NULL;
                        delete altResult;
                        altResult = NULL;
                    }
                }
                if(supportSrsID.size()){
                    auto iit = supportSrsID.begin();
                    int32_t ixmx = iit->first;
                    // first run
                    for(; iit != supportSrsID.end(); ++iit){
                        if(svs[iit->first].mSRSupport > svs[ixmx].mSRSupport){
                            ixmx = iit->first;
                        }
                    }
                    // second run
                    iit = supportSrsID.begin();
                    for(; iit != supportSrsID.end(); ++iit){
                        if(svs[iit->first].mPESupport > svs[ixmx].mPESupport){
                            ixmx = iit->first;
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
        if(!(b->core.flag & BAM_FPAIRED)) continue;
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
            if((sst != 2 || outerISize < mOpt->libInfo->mMinNormalISize || outerISize > mOpt->libInfo->mMaxNormalISize) || (b->core.tid != b->core.mtid)){
                // Get SV type
                if(svt == -1) continue;
                // Spanning a breakpoint?
                bool spanvalid = false;
                int32_t pbegin = b->core.pos;
                int32_t pend = std::min((int32_t)(b->core.pos + mOpt->libInfo->mMaxNormalISize), (int32_t)h->target_len[refIdx]);
                if(b->core.flag & BAM_FREVERSE){
                    pbegin = std::max((BIGD_TYPE)0, b->core.pos + b->core.l_qseq - mOpt->libInfo->mMaxNormalISize);
                    pend = std::min((int32_t)(b->core.pos + b->core.l_qseq), (int32_t)h->target_len[refIdx]);
                }
                for(int32_t i = pbegin; i < pend; ++i){
                    if(spanBp.find(i) != spanBp.end()){
                        spanvalid = true;
                        break;
                    }
                }
                if(spanvalid){
                    std::map<int32_t, bool> supportSpnID;
                    if(lspap != pbegin){
                        itspna = std::lower_bound(spPts[refIdx].begin(), spPts[refIdx].end(), SpanPoint(pbegin));
                        lspap = pbegin;
                        ltap = itspna;
                    }else{
                        itspna = ltap;
                    }
                    for(; itspna != spPts[refIdx].end() && pend >= itspna->mBpPos; ++itspna){
                        if(svt == itspna->mSVT && svs[itspna->mID].mChr1 == b->core.tid && svs[itspna->mID].mChr2 == b->core.mtid){
                            // valid dp in pe-mode
                            bool validDPE = true;
                            if(itspna->mIsSVEnd){
                                if(std::abs(b->core.pos - svs[itspna->mID].mSVEnd) > mOpt->libInfo->mMaxNormalISize){
                                    validDPE = false;
                                }
                                if(std::abs(b->core.mpos - svs[itspna->mID].mSVStart) > mOpt->libInfo->mMaxNormalISize){
                                    validDPE = false;
                                }
                            }else{
                                if(std::abs(b->core.pos - svs[itspna->mID].mSVStart) > mOpt->libInfo->mMaxNormalISize){
                                    validDPE = false;
                                }
                                if(std::abs(b->core.mpos - svs[itspna->mID].mSVEnd) > mOpt->libInfo->mMaxNormalISize){
                                    validDPE = false;
                                }
                            }
                            if(validDPE){
                                supportSpnID[itspna->mID] = itspna->mIsSVEnd;
                            }
                        }
                    }
                    if(supportSpnID.size()){
                        auto iit = supportSpnID.begin();
                        int32_t ixmx = iit->first;
                        // first run
                        for(; iit != supportSpnID.end(); ++iit){
                            if(svs[iit->first].mSRSupport > svs[ixmx].mSRSupport){
                                ixmx = iit->first;
                            }
                        }
                        // second run
                        iit = supportSpnID.begin();
                        for(; iit != supportSpnID.end(); ++iit){
                            if(svs[iit->first].mPESupport > svs[ixmx].mPESupport){
                                ixmx = iit->first;
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
