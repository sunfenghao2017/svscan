#include "svutil.h"
#include "svscanner.h"
#include "ThreadPool.h"

void SVScanner::scanDPandSROne(int32_t tid, JunctionMap* jctMap, DPBamRecordSet* dprSet){
    // Open file handles
    samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
    hts_idx_t* idx = sam_index_load(fp, mOpt->bamfile.c_str());
    hts_set_fai_filename(fp, mOpt->genome.c_str());
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    const uint16_t BAM_SRSKIP_MASK = (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY);
    // Iterate bam contig by contig
    if(mValidRegs[tid].empty()) return; // Skip invalid contig
    uint64_t mapped = 0;
    uint64_t unmapped = 0;
    hts_idx_get_stat(idx, tid, &mapped, &unmapped);
    if(!mapped) return; // Skip contig without any mapped reads
    // Iterate all read alignments on this contig and valid regions
    util::loginfo("Contig: " + std::string(h->target_name[tid]) + " starts SR and DP scanning", mOpt->logMtx);
    for(auto regit = mValidRegs[tid].begin(); regit != mValidRegs[tid].end(); ++regit){
        hts_itr_t* itr = sam_itr_queryi(idx, tid, regit->first, regit->second);
        while(sam_itr_next(fp, itr, b) >= 0){
            if(b->core.flag & BAM_SRSKIP_MASK) continue;// skip invalid reads
            if(b->core.qual < mOpt->filterOpt->minMapQual || b->core.tid < 0) continue;// skip quality poor read
            // SR parsing
            jctMap->insertJunction(b, h);
            // DP parsing
            if(mOpt->libInfo->mMedian == 0) continue; // skip SE library
            if(b->core.flag & BAM_FMUNMAP) continue;// skip invalid reads
            if(mValidRegs[b->core.mtid].empty()) continue;// skip invalid regions
            if(b->core.tid != b->core.mtid && b->core.qual < mOpt->filterOpt->mMinTraQual) continue;// skip quality poor read
            int32_t svt = DPBamRecord::getSVType(b, mOpt);// get sv type
            if(svt == -1) continue; // Skip PE which does not support any SV
            if(mOpt->SVTSet.find(svt) == mOpt->SVTSet.end()) continue;// Skip SV type which does not needed to called
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
    util::loginfo("Start scanning bam for SRs and DPs");
    // Parallel processing all contigs
    samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    int32_t ntargets = h->n_targets;
    // Preset SR and DP result storing objects
    std::vector<JunctionMap*> jct(ntargets, NULL);
    std::vector<DPBamRecordSet*> dps(ntargets, NULL);
    for(int32_t i = 0; i < ntargets; ++i){
        jct[i] = new JunctionMap(mOpt);
        dps[i] = new DPBamRecordSet(mOpt);
    }
    // Construct thread pool
    ThreadPool::ThreadPool pool(ntargets);
    std::vector<std::future<void>> scanret(ntargets);
    for(int32_t i = 0; i < ntargets; ++i){
        scanret[i] = pool.enqueue(&SVScanner::scanDPandSROne, this, i, jct[i], dps[i]);
    }
    for(auto& e: scanret) e.get();
    // Merge and clean
    JunctionMap* jctMap = JunctionMap::merge(jct, mOpt);
    DPBamRecordSet* dprSet = DPBamRecordSet::merge(dps, mOpt);
    for(int32_t i = 0; i < ntargets; ++i){
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
    util::loginfo("Finish scanning bam for SRs and DPs");
    SRBamRecordSet srs(mOpt, jctMap);
    // Update all valid Ref ID s of SR
    for(auto& e: srs.mSRs){
        for(auto& f: e){
            mOpt->svRefID.insert(f.mChr1);
            mOpt->svRefID.insert(f.mChr2);
        }
    }
    util::loginfo("Start clustering SRs");
    srs.cluster(mSRSVs);
    util::loginfo("Finish clustering SRs");
    util::loginfo("Start assembling SRs and refining breakpoints");
    srs.assembleSplitReads(mSRSVs);
    util::loginfo("Finish assembling SRs and refining breakpoints");
    util::loginfo("Found SRSV Candidates: " + std::to_string(mSRSVs.size()));
    // Process all DPs
    util::loginfo("Start clustering DPs");
    dprSet->cluster(mDPSVs);
    util::loginfo("Finish clustering DPs");
    util::loginfo("Found DPSV Candidates: " + std::to_string(mDPSVs.size()));
    // Merge SR and DP SVs
    util::loginfo("Start merging SVs from SRs and DPs");
    SVSet mergedSVs;
    mergeAndSortSVSet(mSRSVs, mDPSVs, mergedSVs);
    util::loginfo("Finish merging SVs from SRs and DPs");
    util::loginfo("Start fetching reference of SV supported by DP only");
    getDPSVRef(mergedSVs, mOpt);
    std::sort(mergedSVs.begin(), mergedSVs.end());
    util::loginfo("Finish fetching reference of SV supported by DP only");
    // Get Allele info of SVs
    for(uint32_t i = 0; i < mergedSVs.size(); ++i){
        mergedSVs[i].addAlleles();
        mergedSVs[i].mID = i;
    }
    // open bamout for write
    mOpt->fbamout = sam_open(mOpt->bamout.c_str(), "w");
    assert(sam_hdr_write(mOpt->fbamout, h) >= 0);
    // Annotate junction reads and spaning coverage
    util::loginfo("Start annotating SV coverage");
    Annotator* covAnn = new Annotator(mOpt);
    Stats* covStat = covAnn->covAnnotate(mergedSVs);
    util::loginfo("Finish annotating SV coverage");
    sam_close(mOpt->fbamout);
    GeneInfoList gl;
    util::loginfo("Start annotating SV gene information");
    covAnn->geneAnnotate(mergedSVs, gl);
    util::loginfo("Finish annotating SV gene information");
    util::loginfo("Start writing SVs to TSV file");
    covStat->reportTSV(mergedSVs, gl);
    util::loginfo("Finish writing SVs to TSV file");
    util::loginfo("Start writing SVs to BCF file");
    covStat->reportBCF(mergedSVs);
    util::loginfo("Finish writing SVs to BCF file");
    delete covAnn;
    delete covStat;
    sam_close(fp);
    bam_hdr_destroy(h);
}
