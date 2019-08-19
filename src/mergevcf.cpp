#include "mergevcf.h"

void VCFMerger::getContigMap(ContigMap& cmap){
    int32_t numseq = 0;
    for(uint32_t c = 0; c < infiles.size(); ++c){
        htsFile* fp = bcf_open(infiles[c].c_str(), "r");
        bcf_hdr_t* hdr = bcf_hdr_read(fp);
        int32_t nseqs = -1;
        const char** seqnames= bcf_hdr_seqnames(hdr, &nseqs);
        for(int32_t i = 0; i < nseqs; ++i){
            if(cmap.find(seqnames[i]) == cmap.end()) cmap[seqnames[i]] = numseq++;
        }
        if(seqnames) free(seqnames);
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
    }
}
void VCFMerger::getIntervals(ContigIntervals& ci, ContigMap& cmap, int32_t insvt){
    for(uint32_t c = 0; c < infiles.size(); ++c){
        htsFile* ifp = bcf_open(infiles[c].c_str(), "r");
        bcf_hdr_t* hdr = bcf_hdr_read(ifp);
        bcf1_t* rec = bcf_init1();
        int32_t nsvend = 0;
        int32_t* svend = NULL;
        int32_t ninslen = 0;
        int32_t* inslen = NULL;
        int32_t npe = 0;
        int32_t* pe = NULL;
        int32_t nsr = 0;
        int32_t* sr = NULL;
        int32_t npemapq = 0;
        int32_t* pemapq = NULL;
        int32_t nsrmapq = 0;
        int32_t* srmapq = NULL;
        int32_t nct = 0;
        char* ct = NULL;
        int32_t nsvt = 0;
        char* svt = NULL;
        while(bcf_read(ifp, hdr, rec) >= 0){
            bcf_unpack(rec, BCF_UN_INFO);
            // skip record didn't PASS
            if(filterForPass && bcf_has_filter(hdr, rec, const_cast<char*>("PASS")) != 1) continue;
            // convert svt back to integer
            int32_t recsvt = -1;
            if(bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt) > 0 && bcf_get_info_string(hdr, rec, "CT", &ct, &nct) > 0){
                recsvt = svutil::str2svt(ct, svt);
            }
            // skip other sv types
            if(recsvt != insvt) continue;
            // unique tid
            int32_t tid = cmap[bcf_hdr_id2name(hdr, rec->rid)];
            // sv start and end
            int32_t svStart = rec->pos;
            int32_t svEnd = svStart + 1;
            if(bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) svEnd = *svend;
            // sv size/check
            if(std::strcmp(svt, "BND") != 0){
                if(std::strcmp(svt, "INS") == 0){
                    if(bcf_get_info_int32(hdr, rec, "INSLEN", &inslen, &ninslen) > 0){
                        if(*inslen < minsize || *inslen > maxsize) continue;// skip invalid sv size
                    }
                }else{
                    if(svEnd - svStart < minsize || svEnd - svStart > maxsize) continue;// skip invalid sv size
                }
            }
            // skip record which is not PRECISE if needed
            bool precise = false;
            if(bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) == 1) precise = true;
            if(filterForPrecise && (!precise)) continue;
            // variant allele frequency filter
            if(vaf > 0 || coverage > 0){
                float maxvaf = 0;
                int32_t maxcov = 0;
                bcf_unpack(rec, BCF_UN_ALL);
                int32_t ndv = 0;
                int32_t* dv = NULL;
                int32_t ndr = 0;
                int32_t* dr = 0;
                int32_t nrv = 0;
                int32_t* rv = NULL;
                int32_t nrr = 0;
                int32_t* rr = 0;
                int ngt = 0;
                int32_t* gt = NULL;
                bcf_get_format_int32(hdr, rec, "DV", &dv, &ndv);
                bcf_get_format_int32(hdr, rec, "DR", &dr, &ndr);
                bcf_get_format_int32(hdr, rec, "RV", &rv, &nrv);
                bcf_get_format_int32(hdr, rec, "RR", &rr, &nrr);
                bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
                for(int32_t i = 0; i < bcf_hdr_nsamples(hdr); ++i){
                    int32_t supportSum = 0;
                    if(precise) supportSum = rr[i] + rv[i];
                    else supportSum = dr[i] + dv[i];
                    if(supportSum > 0){
                        double vaf = 0;
                        if(precise) vaf = (double)rv[i]/(double)supportSum;
                        else vaf = (double)dv[i]/(double)supportSum;
                        if(vaf > maxvaf) maxvaf = vaf;
                        if(supportSum > maxcov) maxcov = supportSum;
                    }
                }
                if(dv) free(dv);
                if(dr) free(dr);
                if(rv) free(rv);
                if(rr) free(rr);
                if(gt) free(gt);
                if(maxvaf < vaf || maxcov < coverage) continue;
            }
            // get sv support and mapq
            int32_t peSupport = 0, srSupport = 0;
            int32_t peMapQual = 0, srMapQual = 0;
            if(bcf_get_info_int32(hdr, rec, "PE", &pe, &npe) > 0) peSupport = *pe;
            if(bcf_get_info_int32(hdr, rec, "SR", &sr, &nsr) > 0) srSupport = *sr;
            if(bcf_get_info_int32(hdr, rec, "PEMAPQ", &pemapq, &npemapq) > 0) peMapQual = *pemapq;
            if(bcf_get_info_int32(hdr, rec, "SRMAPQ", &srmapq, &nsrmapq) > 0) srMapQual = *srmapq;
            // quality score for the sv
            int32_t svScore = 3 * srSupport * srMapQual + peSupport * peMapQual;
            if(bcf_hdr_id2int(hdr, BCF_DT_ID, "SCORE") >= 0){
                int32_t nscore = 0;
                if(bcf_hdr_id2type(hdr, BCF_HL_INFO, bcf_hdr_id2int(hdr, BCF_DT_ID, "SCORE")) == BCF_HT_INT){
                    int32_t* score = NULL;
                    bcf_get_info_int32(hdr, rec, "SCORE", &score, &nscore);
                    svScore = *score;
                    free(score);
                }else if(bcf_hdr_id2type(hdr, BCF_HL_INFO, bcf_hdr_id2int(hdr, BCF_DT_ID, "SCORE")) == BCF_HT_REAL){
                    float* score = NULL;
                    bcf_get_info_float(hdr, rec, "SCORE", & score, &nscore);
                    svScore = 10000 * *score;
                    free(score);
                }
            }
            // store the interval
            ci[tid].push_back(Interval(svStart, svEnd, svScore));
        }
        if(svend) free(svend);
        if(inslen) free(inslen);
        if(pe) free(pe);
        if(sr) free(sr);
        if(pemapq) free(pemapq);
        if(srmapq) free(srmapq);
        if(svt) free(svt);
        bcf_hdr_destroy(hdr);
        bcf_destroy(rec);
        bcf_close(ifp);
    }
    // sort Interval by sv start and end pos
    for(uint32_t i = 0; i < ci.size(); ++i) std::sort(ci[i].begin(), ci[i].end());
}

void VCFMerger::mergeIntervals(ContigIntervals& ici, ContigIntervals& oci, int32_t svt){
    int32_t tid = 0;
    // merge
    for(auto it = ici.begin(); it != ici.end(); ++it, ++tid){
        std::vector<bool> keepMarker(it->size(), true);
        auto ik = keepMarker.begin();
        for(auto is = it->begin(); is != it->end(); ++is, ++ik){
            auto isNext = is;
            auto ikNext = ik;
            ++isNext;
            ++ikNext;
            for(; isNext != it->end(); ++isNext, ++ikNext){
                if(isNext->mStart - is->mStart > bpoffset) break;
                else{
                    if(std::abs(isNext->mEnd - is->mEnd) < bpoffset){
                        if(svt >= 5 || is->overlap(*isNext) >= recoverlap){
                            if(is->mScore < isNext->mScore) *ik = false;
                            else if(isNext->mScore < is->mScore) *ikNext = false;
                            else{
                                if(is->mStart < isNext->mStart) *ikNext = false;
                                else if(is->mEnd < isNext->mEnd) *ikNext = false;
                                else *ik = false;
                            }
                        }
                    }
                }
            }
            if(*ik) oci[tid].push_back(Interval(is->mScore, is->mEnd, is->mScore));
        }
    }
    // sort
    for(int32_t i = 0; i < tid; ++i)  std::sort(oci[i].begin(), oci[i].end());
}
