#include "svutil.h"
#include "svscanner.h"
#include "ThreadPool.h"

void SVScanner::scanDPandSROne(int32_t tid, JunctionMap* jctMap, DPBamRecordSet* dprSet){
    if(mScanRegs[tid].empty()) return; // Skip invalid contig
    // Open file handles
    samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
    hts_idx_t* idx = sam_index_load(fp, mOpt->bamfile.c_str());
    hts_set_fai_filename(fp, mOpt->alnref.c_str());
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    const uint16_t BAM_SRSKIP_MASK = (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY);
    // Iterate bam contig by contig
    uint64_t mapped = 0;
    uint64_t unmapped = 0;
    hts_idx_get_stat(idx, tid, &mapped, &unmapped);
    if(!mapped){
        sam_close(fp);
        hts_idx_destroy(idx);
        bam_hdr_destroy(h);
        bam_destroy1(b);
        return; // Skip contig without any mapped reads
    }
    // Iterate all read alignments on this contig and valid regions
    util::loginfo("Contig: " + std::string(h->target_name[tid]) + " starts SR and DP scanning", mOpt->logMtx);
    for(auto regit = mScanRegs[tid].begin(); regit != mScanRegs[tid].end(); ++regit){
        hts_itr_t* itr = sam_itr_queryi(idx, tid, regit->first, regit->second);
        while(sam_itr_next(fp, itr, b) >= 0){
            if(b->core.flag & BAM_SRSKIP_MASK) continue;// skip invalid reads
            if(b->core.qual < mOpt->filterOpt->minMapQual || b->core.tid < 0) continue;// skip quality poor read
            if(!inValidReg(b, h)) continue; // skip reads which do not overlap with creg, neither does its mate
            int lhc = 0, rhc = 0;
            std::pair<int, int> sclens = bamutil::getClipLength(b, lhc, rhc);
            if(lhc || rhc) continue; // skip reads wich any hard clips
            if(sclens.first && sclens.second) continue; // skip reads with heading and leading clips
            if(sclens.first + sclens.second) jctMap->insertJunction(b, h); // only one softclip read can be SR candidates
            if(mOpt->libInfo->mMedian == 0) continue; // skip SE library from DP collecting
            if(b->core.flag & BAM_FMUNMAP) continue;// skip invalid reads
            if(mScanRegs[b->core.mtid].empty()) continue;// skip invalid regions
            if(b->core.tid != b->core.mtid && b->core.qual < mOpt->filterOpt->mMinTraQual) continue;// skip quality poor read
            int32_t svt = DPBamRecord::getSVType(b, mOpt);// get sv type
            if(svt == -1) continue; // Skip PE which does not support any SV
            if(mOpt->SVTSet.find(svt) == mOpt->SVTSet.end()) continue;// Skip SV type which does not needed to called
            uint8_t* data = bam_aux_get(b, "MC");
            if(data){
                char* c = bam_aux2Z(data);
                int32_t leadingSC = 0, tailingSC = 0;
                while(*c && *c != '*'){
                    int32_t num = 0;
                    if(std::isdigit((int)*c)){
                        num = std::strtol(c, &c, 10);
                    }else{
                        num = 1;
                    }
                    switch(*c){
                        case 'S':
                            if(!leadingSC) leadingSC = num;
                            else tailingSC = num;
                            break;
                        default:
                            break;
                    }
                    ++c;
                }
                if(leadingSC && tailingSC) continue; // skip mate with both leading/tailing clips
            }
            if(b->core.tid > b->core.mtid || (b->core.tid == b->core.mtid && b->core.pos > b->core.mpos)){// second read in pair
                dprSet->insertDP(b, svt);
            }
        }
        hts_itr_destroy(itr);
    }
    util::loginfo("Contig: " + std::string(h->target_name[tid]) + " finished SR and DP scanning", mOpt->logMtx);
    sam_close(fp);
    hts_idx_destroy(idx);
    bam_hdr_destroy(h);
    bam_destroy1(b);
}

void SVScanner::scanDPandSR(){
    util::loginfo("Beg scanning bam for SRs and DPs");
    // Preset SR and DP result storing objects
    std::vector<JunctionMap*> jct(mOpt->contigNum, NULL);
    std::vector<DPBamRecordSet*> dps(mOpt->contigNum, NULL);
    for(int32_t i = 0; i < mOpt->contigNum; ++i){
        jct[i] = new JunctionMap(mOpt);
        dps[i] = new DPBamRecordSet(mOpt);
    }
    // Stat reads count in each contig
    samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
    hts_idx_t* idx = sam_index_load(fp, mOpt->bamfile.c_str());
    std::vector<RegItemCnt> ctgRdStat;
    for(int32_t i = 0; i < mOpt->contigNum; ++i){
        uint64_t unmapped  = 0, mapped = 0;
        if(hts_idx_get_stat(idx, i, &mapped, &unmapped) >= 0){
            if(mapped > 0){
                RegItemCnt crc;
                crc.mTid = i;
                crc.mCount = mapped;
                ctgRdStat.push_back(crc);
            }
        }
    }
    std::sort(ctgRdStat.begin(), ctgRdStat.end());
    sam_close(fp);
    hts_idx_destroy(idx);
    // Parallel processing each contig
    std::vector<std::future<void>> scanret(ctgRdStat.size());
    for(uint32_t i = 0; i < ctgRdStat.size(); ++i){
        int32_t refidx = ctgRdStat[i].mTid;
        scanret[i] = mOpt->pool->enqueue(&SVScanner::scanDPandSROne, this, refidx, jct[refidx], dps[refidx]);
    }
    for(auto& e: scanret) e.get();
    // Merge and clean
    JunctionMap* jctMap = JunctionMap::merge(jct, mOpt);
    DPBamRecordSet* dprSet = DPBamRecordSet::merge(dps, mOpt);
    for(int32_t i = 0; i < mOpt->contigNum; ++i){
        delete jct[i];
        jct[i] = NULL;
        delete dps[i];
        dps[i] = NULL;
    }
    // Update all valid Ref IDs of DP
    for(auto& e: dprSet->mDPs){
        for(auto& f: e){
            mOpt->svRefID.insert(f.mCurTid);
            mOpt->svRefID.insert(f.mMateTid);
        }
    }
    // Process all SRs
    util::loginfo("End scanning bam for SRs and DPs");
    SRBamRecordSet srs(mOpt, jctMap);
    util::loginfo("Beg clustering SRs");
    // Update all valid Ref IDs of SR
    for(auto& e: srs.mSRs){
        for(auto& f: e){
            mOpt->svRefID.insert(f.mChr1);
            mOpt->svRefID.insert(f.mChr2);
        }
    }
    if(mOpt->debug & DEBUG_FCALL){
        std::cout << "debug_svRefID:";
        for(auto& srfid: mOpt->svRefID){
            std::cout << "\t" << srfid;
        }
        std::cout << std::endl;
    }
    srs.cluster(mSRSVs);
    util::loginfo("End clustering SRs");
    if(mOpt->debug & DEBUG_FCALL){
        std::cout << "debug_SRSV_Cluster_result: " << std::endl;
        std::cout << srs << std::endl;
    }
    util::loginfo("Beg assembling SRs and refining breakpoints");
    srs.assembleSplitReads(mSRSVs);
    util::loginfo("End assembling SRs and refining breakpoints");
    util::loginfo("Found SRSV Candidates: " + std::to_string(mSRSVs.size()));
    // Process all DPs
    util::loginfo("Beg clustering DPs");
    dprSet->cluster(mDPSVs);
    util::loginfo("End clustering DPs");
    if(mOpt->debug & DEBUG_FCALL){
        std::cout << "debug_DP_SV_Cluster_result: " << std::endl;
        std::cout << dprSet << std::endl;
    }
    util::loginfo("Found DPSV Candidates: " + std::to_string(mDPSVs.size()));
    // sort all svs
    std::sort(mDPSVs.begin(), mDPSVs.end());
    for(uint32_t i = 0; i < mDPSVs.size(); ++i){
        mDPSVs[i].mID = i;
    }
    std::sort(mSRSVs.begin(), mSRSVs.end());
    for(uint32_t i = 0; i < mSRSVs.size(); ++i){
        mSRSVs[i].mID = i;
    }
    if(mOpt->debug & DEBUG_FCALL){
        std::cout << "debug_DP_SVS_found: " << std::endl;
        std::cout << mDPSVs << std::endl;
        std::cout << "debug_SR_SVS_found: " << std::endl;
        std::cout << mSRSVs << std::endl;
    }
    // Merge SR and DP SVs
    util::loginfo("Beg merging SVs from SRs and DPs");
    SVSet mergedSVs;
    mergeAndSortSVSet(mSRSVs, mDPSVs, mergedSVs, mOpt);
    if(mOpt->debug & DEBUG_FCALL){
        std::cout << "debug_ALL_SVS_found: " << std::endl;
        std::cout << mergedSVs << std::endl;
    }
    util::loginfo("End merging SVs from SRs and DPs, all SV got: " + std::to_string(mergedSVs.size()));
    util::loginfo("Beg fetching reference of SV supported by DP only");
    getDPSVRef(mergedSVs, mOpt);
    std::sort(mergedSVs.begin(), mergedSVs.end());
    util::loginfo("End fetching reference of SV supported by DP only");
    // Get Allele info of SVs and updatev SVsupporting contigs as well as svsize
    mOpt->svRefID.clear();
    for(uint32_t i = 0; i < mergedSVs.size(); ++i){
        mergedSVs[i].addAlleles();
        mergedSVs[i].mID = i;
        if(mergedSVs[i].mSVT >= 5){
            mergedSVs[i].mSize = -1;
        }else if(mergedSVs[i].mSVT == 4){
            mergedSVs[i].mSize = mergedSVs[i].mConsensus.size();
        }else{
            mergedSVs[i].mSize = mergedSVs[i].mSVEnd - mergedSVs[i].mSVStart;
        }
        mOpt->svRefID.insert(mergedSVs[i].mChr1);
        mOpt->svRefID.insert(mergedSVs[i].mChr2);
    }
    // open bamout for write
    if(!mOpt->bamout.empty()){
        samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
        bam_hdr_t* h = sam_hdr_read(fp);
        mOpt->fbamout = sam_open(mOpt->bamout.c_str(), "w");
        assert(sam_hdr_write(mOpt->fbamout, h) >= 0);
        sam_close(fp);
        bam_hdr_destroy(h);
    }
    // Annotate junction reads and spaning coverage
    util::loginfo("Beg annotating SV coverage");
    Annotator* covAnn = new Annotator(mOpt);
    Stats* covStat = covAnn->covAnnotate(mergedSVs);
    util::loginfo("End annotating SV coverage");
    if(!mOpt->bamout.empty()) sam_close(mOpt->fbamout);
    GeneInfoList gl;
    util::loginfo("Beg annotating SV gene information");
    if(mOpt->rnamode){
        covAnn->geneAnnoRNA(mergedSVs, gl);
    }else{
        covAnn->geneAnnoDNA(mergedSVs, gl);
    }
    covStat->maskFuseRec(mergedSVs, gl);
    util::loginfo("End annotating SV gene information");
    util::loginfo("Beg writing SVs to TSV file");
    covStat->reportSVTSV(mergedSVs, gl);
    util::loginfo("End writing SVs to TSV file");
    util::loginfo("Beg writing Fusions to TSV file");
    covStat->reportFusionTSV(mergedSVs, gl);
    util::loginfo("End writing Fusions to TSV file");
    if(!mOpt->bcfOut.empty()){
        util::loginfo("Beg writing SVs to BCF file");
        covStat->reportSVBCF(mergedSVs);
        util::loginfo("End writing SVs to BCF file");
    }
    if(!mOpt->bam2tb.empty() && !mOpt->bamout.empty()){
        util::loginfo("Beg writing fusion supporting bam records to excel file");
        BamToTable btt;
        btt.svbam = mOpt->bamout;
        btt.fstsv = mOpt->fuseOpt->mOutFile;
        btt.sstsv = mOpt->fuseOpt->mSupFile;
        btt.bamtb = mOpt->bam2tb;
        btt.b2t();
        util::loginfo("End writing fusion supporting bam records to excel file");
    }
    delete covAnn;
    delete covStat;
}
