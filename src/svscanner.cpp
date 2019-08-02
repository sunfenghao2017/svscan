#include "svutil.h"
#include "svscanner.h"

void SVScanner::scanDPandSR(){
    // Open file handles
    samFile* fp = sam_open(mOpt->bamfile.c_str(), "r");
    hts_idx_t* idx = sam_index_load(fp, mOpt->bamfile.c_str());
    hts_set_fai_filename(fp, mOpt->genome.c_str());
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    // Split-read records
    JunctionMap* jctMap = new JunctionMap(mOpt);
    const uint16_t BAM_SRSKIP_MASK = (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY);
    // DP-read records
    DPBamRecordSet* dprSet = new DPBamRecordSet(mOpt);
    const uint16_t BAM_DPSKIP_MASK = (BAM_FSUPPLEMENTARY | BAM_FMUNMAP);
    // Inter-chromosomal mate map and alignment length
    std::unordered_map<size_t, std::pair<uint8_t, int32_t>> matetra;
    std::unordered_map<size_t, std::pair<uint8_t, int32_t>> matemap;
    util::loginfo("Start scanning bam for SRs and DPs");
    // Iterate bam contig by contig
    for(int32_t refIndex = 0; refIndex < h->n_targets; ++refIndex){
        if(mValidRegs[refIndex].empty()) continue; // Skip invalid contig
        uint64_t mapped = 0;
        uint64_t unmapped = 0;
        hts_idx_get_stat(idx, refIndex, &mapped, &unmapped);
        if(!mapped) continue; // Skip contig without any mapped reads
        // Iterate all read alignments on this contig and valid regions
        for(auto regit = mValidRegs[refIndex].begin(); regit != mValidRegs[refIndex].end(); ++regit){
            hts_itr_t* itr = sam_itr_queryi(idx, refIndex, regit->first, regit->second);
            int32_t lastAlignedPos = 0;
            std::set<size_t> lastAlignedPosReads;
            while(sam_itr_next(fp, itr, b) >= 0){
                if(b->core.flag & BAM_SRSKIP_MASK) continue;// skip invalid reads
                if(b->core.qual < mOpt->filterOpt->minMapQual || b->core.tid < 0) continue;// skip quality poor read
                // Try to parse and insert an SR bam record
                if(jctMap->insertJunction(b)) mOpt->svRefID.insert(b->core.tid);
                // DP parsing
                if(mOpt->libInfo->mMedian == 0) continue; // skip SE library
                if(b->core.flag & BAM_DPSKIP_MASK) continue;// skip invalid reads
                if(mValidRegs[b->core.mtid].empty()) continue;// skip invalid regions
                if(b->core.tid != b->core.mtid && b->core.qual < mOpt->filterOpt->mMinTraQual) continue;// skip quality poor read
                int32_t svt = DPBamRecord::getSVType(b, mOpt);// get sv type
                if(svt == -1) continue; // Skip PE which does not support any SV
                if(mOpt->SVTSet.find(svt) == mOpt->SVTSet.end()) continue;// Skip SV type which does not needed to called
                if(b->core.pos > lastAlignedPos){// clear records aligned at the same position
                    lastAlignedPosReads.clear();
                    lastAlignedPos = b->core.pos;
                }
                if(Stats::firstInPair(b, lastAlignedPosReads)){// First in pair
                    size_t hv = svutil::hashPairCurr(b);
                    lastAlignedPosReads.insert(hv);
                    if(svt >= 5) matetra[hv] = std::make_pair(b->core.pos, bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)));
                    else matemap[hv] = std::make_pair(b->core.pos, bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)));
                }else{// Second in pair
                    size_t hv = svutil::hashPairMate(b);
                    int32_t matealn = 0;
                    uint8_t pairQual = 0;
                    if(svt >= 5){// translocation
                        auto mit = matetra.find(hv);
                        if(mit == matetra.end()) continue; // Skip read whose mate discarded
                        if(mit->second.first == 0) continue; // Skip read whose mate is mapped to multiple place
                        pairQual = std::min(mit->second.first, b->core.qual);
                        matealn = mit->second.second;
                        mit->second.first = 0;
                    }else{// sv on same chr
                        auto mit = matemap.find(hv);
                        if(mit == matemap.end()) continue; // Skip read whose mate discarded
                        if(mit->second.first == 0) continue; // Skip read whose mate is mapped to multiple place
                        pairQual = std::min(mit->second.first, b->core.qual);
                        matealn = mit->second.second;
                        mit->second.first = 0;
                    }
                    dprSet->insertDP(b, matealn, pairQual, svt);
                    mOpt->svRefID.insert(b->core.tid);
                    mOpt->svRefID.insert(b->core.mtid);
                    ++mOpt->libInfo->mAbnormalPairs;
                }
            }
            hts_itr_destroy(itr);
        }
        util::loginfo("Contig: " + std::string(h->target_name[refIndex]) + " finished SR and DP scanning");
    }
    bam_destroy1(b);
    // Process all SRs
    util::loginfo("Finish scanning bam for SRs and DPs");
    SRBamRecordSet srs(mOpt, jctMap);
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
    mergeAndSortSVSet(mSRSVs, mDPSVs, 500, 10);
    util::loginfo("Finish merging SVs from SRs and DPs");
    util::loginfo("Start fetching reference of SV supported by DP only");
    getDPSVRef(mDPSVs, mOpt);
    util::loginfo("Finish fetching reference of SV supported by DP only");
    // Get Allele info of SVs
    for(uint32_t i = 0; i < mDPSVs.size(); ++i) mDPSVs[i].addAlleles();
    // Annotate junction reads and spaning coverage
    util::loginfo("Start annotating SVs");
    Annotator* covAnn = new Annotator(mOpt);
    Stats* covStat = covAnn->covAnnotate(mDPSVs);
    util::loginfo("Finish annotating SVs");
    util::loginfo("Start writing SVs to BCF file");
    covStat->reportBCF(mDPSVs);
    util::loginfo("Finish writing SVs to BCF file");
    delete covAnn;
    delete covStat;
}
