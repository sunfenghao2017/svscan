#include "stats.h"

Stats::Stats(Options* opt, int32_t n){
    mOpt = opt;
    init(n);
}

void Stats::init(int n){
    mJctCnts.resize(n, JunctionCount());
    mSpnCnts.resize(n, SpanningCount());
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
            // SC
            ret->mJctCnts[j].mAlth1 += sts[i]->mJctCnts[j].mAlth1;
            ret->mJctCnts[j].mAlth2 += sts[i]->mJctCnts[j].mAlth2;
            ret->mJctCnts[j].mRefh1 += sts[i]->mJctCnts[j].mRefh1;
            ret->mJctCnts[j].mFPIns += sts[i]->mJctCnts[j].mFPIns;
            ret->mJctCnts[j].mAltCnt += sts[i]->mJctCnts[j].mAltCnt;
            ret->mJctCnts[j].mRefCntBeg += sts[i]->mJctCnts[j].mRefCntBeg;
            ret->mJctCnts[j].mRefCntEnd += sts[i]->mJctCnts[j].mRefCntEnd;
            ret->mJctCnts[j].mAltQual.insert(ret->mJctCnts[j].mAltQual.end(), sts[i]->mJctCnts[j].mAltQual.begin(), sts[i]->mJctCnts[j].mAltQual.end());
            ret->mJctCnts[j].mRefQualBeg.insert(ret->mJctCnts[j].mRefQualBeg.end(), sts[i]->mJctCnts[j].mRefQualBeg.begin(), sts[i]->mJctCnts[j].mRefQualBeg.end());
            ret->mJctCnts[j].mRefQualEnd.insert(ret->mJctCnts[j].mRefQualEnd.end(), sts[i]->mJctCnts[j].mRefQualEnd.begin(), sts[i]->mJctCnts[j].mRefQualEnd.end());
            // DP
            ret->mSpnCnts[j].mAlth1 += sts[i]->mSpnCnts[j].mAlth1;
            ret->mSpnCnts[j].mAlth2 += sts[i]->mSpnCnts[j].mAlth2;
            ret->mSpnCnts[j].mRefh1 += sts[i]->mSpnCnts[j].mRefh1;
            ret->mSpnCnts[j].mAltCnt += sts[i]->mSpnCnts[j].mAltCnt;
            ret->mSpnCnts[j].mRefCntBeg += sts[i]->mSpnCnts[j].mRefCntBeg;
            ret->mSpnCnts[j].mRefCntEnd += sts[i]->mSpnCnts[j].mRefCntEnd;
            ret->mSpnCnts[j].mAltQual.insert(ret->mSpnCnts[j].mAltQual.end(), sts[i]->mSpnCnts[j].mAltQual.begin(), sts[i]->mSpnCnts[j].mAltQual.end());
            ret->mSpnCnts[j].mRefQualBeg.insert(ret->mSpnCnts[j].mRefQualBeg.end(), sts[i]->mSpnCnts[j].mRefQualBeg.begin(), sts[i]->mSpnCnts[j].mRefQualBeg.end());
            ret->mSpnCnts[j].mRefQualEnd.insert(ret->mSpnCnts[j].mRefQualEnd.end(), sts[i]->mSpnCnts[j].mRefQualEnd.begin(), sts[i]->mSpnCnts[j].mRefQualEnd.end());
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
                         std::max(0, b->core.pos - mOpt->libInfo->mMaxNormalISize), 
                         std::min(b->core.pos + mOpt->libInfo->mMaxNormalISize, (int32_t)h->target_len[b->core.tid]))){
           continue;
        }
        // Count aligned basepair (small InDels)
        int32_t leadingSC = 0;
        int32_t tailingSC = 0;
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
            }
        }
        if(leadingSC && tailingSC) continue; // skip reads with both leading and tailing softclips
        // Check read length for junction annotation
        if(b->core.l_qseq > 2 * mOpt->filterOpt->mMinFlankSize){
            bool bpvalid = false;
            int32_t rbegin = std::max(0, b->core.pos - leadingSC);
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
                        if(!(leadingSC + tailingSC)){// REF Type, no realignment needed
                            if(b->core.qual >= mOpt->filterOpt->mMinGenoQual){
                                if(itbp->mIsSVEnd){
                                    mJctCnts[itbp->mID].mRefCntEnd += 1;
                                    if(mOpt->writebcf) mJctCnts[itbp->mID].mRefQualEnd.push_back(b->core.qual);
                                }else{
                                    mJctCnts[itbp->mID].mRefCntBeg += 1;
                                    if(mOpt->writebcf) mJctCnts[itbp->mID].mRefQualBeg.push_back(b->core.qual);
                                }
                                uint8_t* hpptr = bam_aux_get(b, "HP");
                                if(hpptr){
                                    mOpt->libInfo->mIsHaploTagged = true;
                                    int hapv = bam_aux2i(hpptr);
                                    if(hapv == 1) ++mJctCnts[itbp->mID].mRefh1;
                                    else ++mJctCnts[itbp->mID].mRefh2; 
                                }
                            }
                            continue;
                        }
                        // possible ALT
                        std::string consProbe = itbp->mIsSVEnd ? svs[itbp->mID].mProbeEndC : svs[itbp->mID].mProbeBegC;
                        std::string refProbe = itbp->mIsSVEnd ? svs[itbp->mID].mProbeEndR : svs[itbp->mID].mProbeBegR;
                        if(readOri.empty()) readOri = bamutil::getSeq(b); // then fetch read to do realign
                        readSeq = readOri;
                        SRBamRecord::adjustOrientation(readSeq, itbp->mIsSVEnd, itbp->mSVT);
                        // Compute alignment to alternative haplotype
                        Aligner* altAligner = new Aligner(consProbe, readSeq, &alnCfg);
                        Matrix2D<char>* altResult = new Matrix2D<char>();
                        int alnScore = altAligner->needle(altResult);
                        int matchThreshold = mOpt->filterOpt->mFlankQuality * consProbe.size() * alnCfg.mMatch + (1 - mOpt->filterOpt->mFlankQuality) * consProbe.size() * alnCfg.mMisMatch;
                        double scoreAlt = (double)alnScore / (double)matchThreshold;
                        // Compute alignment to reference haplotype
                        Aligner* refAligner = new Aligner(refProbe, readOri, &alnCfg);
                        Matrix2D<char>* refResult = new Matrix2D<char>();
                        alnScore = refAligner->needle(refResult);
                        matchThreshold = mOpt->filterOpt->mFlankQuality * refProbe.size() * alnCfg.mMatch + (1 - mOpt->filterOpt->mFlankQuality) * refProbe.size() * alnCfg.mMisMatch;
                        double scoreRef = (double)alnScore / (double)matchThreshold;
                        // Any confident alignment?
                        if(scoreRef > mOpt->filterOpt->mMinSRResScore || scoreAlt > mOpt->filterOpt->mMinSRResScore){
                            if(scoreRef > scoreAlt){
                                if(b->core.qual >= mOpt->filterOpt->mMinGenoQual){
                                    mOpt->logMtx.lock();
                                    if(itbp->mIsSVEnd){
                                        mJctCnts[itbp->mID].mRefCntEnd += 1;
                                        if(mOpt->writebcf) mJctCnts[itbp->mID].mRefQualEnd.push_back(b->core.qual);
                                    }else{
                                        mJctCnts[itbp->mID].mRefCntBeg += 1;
                                        if(mOpt->writebcf) mJctCnts[itbp->mID].mRefQualBeg.push_back(b->core.qual);
                                    }
                                    uint8_t* hpptr = bam_aux_get(b, "HP");
                                    if(hpptr){
                                        mOpt->libInfo->mIsHaploTagged = true;
                                        int hapv = bam_aux2i(hpptr);
                                        if(hapv == 1) ++mJctCnts[itbp->mID].mRefh1;
                                        else ++mJctCnts[itbp->mID].mRefh2;
                                    }
                                    mOpt->logMtx.unlock();
                                }
                            }else{
                                if(itbp->mSVT == 4) supportInsID.push_back(itbp->mID);
                                if(itbp->mSVT != 4) onlySupportIns = false;
                                if(b->core.qual >= mOpt->filterOpt->mMinGenoQual){
                                    mOpt->logMtx.lock();
                                    mJctCnts[itbp->mID].mAltCnt += 1;
                                    if(mOpt->writebcf) mJctCnts[itbp->mID].mAltQual.push_back(b->core.qual);
                                    uint8_t* hpptr = bam_aux_get(b, "HP");
                                    if(hpptr){
                                        mOpt->libInfo->mIsHaploTagged = true;
                                        int hapv = bam_aux2i(hpptr);
                                        if(hapv == 1) ++mJctCnts[itbp->mID].mAlth1;
                                        else ++mJctCnts[itbp->mID].mAlth2;
                                    }
                                    if(mOpt->fbamout){
                                        bam_aux_update_int(b, "ZF", itbp->mID);
                                        assert(sam_write1(mOpt->fbamout, h, b) >= 0);
                                    }
                                    mOpt->logMtx.unlock();
                                }
                            }
                        }
                        delete altAligner;
                        altAligner = NULL;
                        delete altResult;
                        altResult = NULL;
                        delete refAligner;
                        refAligner = NULL;
                        delete refResult;
                        refResult = NULL;
                    }
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
            if((DPBamRecord::getSVType(b) == 2) && outerISize >= mOpt->libInfo->mMinNormalISize &&
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
                            mSpnCnts[itspnr->mID].mRefCntEnd += 1;
                            if(mOpt->writebcf) mSpnCnts[itspnr->mID].mRefQualEnd.push_back(b->core.qual);
                        }else{
                            mSpnCnts[itspnr->mID].mRefCntBeg += 1;
                            if(mOpt->writebcf) mSpnCnts[itspnr->mID].mRefQualBeg.push_back(b->core.qual);
                        }
                        uint8_t* hpptr = bam_aux_get(b, "HP");
                        if(hpptr){
                            mOpt->libInfo->mIsHaploTagged = true;
                            int hap = bam_aux2i(hpptr);
                            if(hap == 1) ++mSpnCnts[itspnr->mID].mRefh1;
                            else ++mSpnCnts[itspnr->mID].mRefh2;
                        }
                        mOpt->logMtx.unlock();
                    }
                }
            }
            // Abnormal spanning coverage
            if(((DPBamRecord::getSVType(b) != 2) || outerISize < mOpt->libInfo->mMinNormalISize || outerISize > mOpt->libInfo->mMaxNormalISize) || (b->core.tid != b->core.mtid)){
                // Get SV type
                int32_t svt =  DPBamRecord::getSVType(b, mOpt);
                if(svt == -1) continue;
                // Spanning a breakpoint?
                bool spanvalid = false;
                int32_t pbegin = b->core.pos;
                int32_t pend = std::min(b->core.pos + mOpt->libInfo->mMaxNormalISize, (int32_t)h->target_len[refIdx]);
                if(b->core.flag & BAM_FREVERSE){
                    pbegin = std::max(0, b->core.pos + b->core.l_qseq - mOpt->libInfo->mMaxNormalISize);
                    pend = std::min(b->core.pos + b->core.l_qseq, (int32_t)h->target_len[refIdx]);
                }
                for(int32_t i = pbegin; i < pend; ++i){
                    if(spanBp.find(i) != spanBp.end()){
                        spanvalid = true;
                        break;
                    }
                }
                if(spanvalid){
                    if(lspap != pbegin){
                        itspna = std::lower_bound(spPts[refIdx].begin(), spPts[refIdx].end(), SpanPoint(pbegin));
                        lspap = pbegin;
                        ltap = itspna;
                    }else{
                        itspna = ltap;
                    }
                    for(; itspna != spPts[refIdx].end() && pend >= itspna->mBpPos; ++itspna){
                        if(svt == itspna->mSVT && svs[itspna->mID].mChr1 == b->core.tid && svs[itspna->mID].mChr2 == b->core.mtid){
                            mOpt->logMtx.lock();
                            mSpnCnts[itspna->mID].mAltCnt += 1;
                            if(mOpt->writebcf) mSpnCnts[itspna->mID].mAltQual.push_back(b->core.qual);
                            uint8_t* hpptr = bam_aux_get(b, "HP");
                            if(hpptr){
                                mOpt->libInfo->mIsHaploTagged = true;
                                int hap = bam_aux2i(hpptr);
                                if(hap == 1) ++mSpnCnts[itspna->mID].mAlth1;
                                else ++mSpnCnts[itspna->mID].mAlth2;
                            }
                            if(mOpt->fbamout){
                                bam_aux_update_int(b, "ZF", itspna->mID);
                                assert(sam_write1(mOpt->fbamout, h, b) >= 0);
                            }
                            mOpt->logMtx.unlock();
                        }
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
