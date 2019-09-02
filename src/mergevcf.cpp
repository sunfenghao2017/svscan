#include "mergevcf.h"

VCFMerger::VCFMerger(){
    softEnv = new Software();
    softEnv->cmp += "version: " + softEnv->version + "\n";
    softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
}

VCFMerger::~VCFMerger(){
    if(softEnv){
        delete softEnv;
    }
}

void VCFMerger::update(int argc, char** argv){
    // update software environment records
    softEnv->cwd = util::cwd();
    for(int i = 0; i < argc; ++i){
        softEnv->cmd.append(argv[i]);
        softEnv->cmd.append(" ");
    }
    // parse all input bcf files to bcf file list
    std::ifstream fr(infilelist);
    std::string fpath;
    while(std::getline(fr, fpath)){
        infiles.push_back(fpath);
    }
    fr.close();
    // check input bcf file loadable
    for(uint32_t c = 0; c < infiles.size(); ++c){
        htsFile* ifp = bcf_open(infiles[c].c_str(), "r");
        if(!ifp) util::errorExit("Failed loading " + infiles[c]);
        bcf_close(ifp);
    }
    // prepare output tempdir
    if(!util::exists(tmpdir)){
        util::makedir(tmpdir);
    }
}

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
            if(*ik) oci[tid].push_back(Interval(is->mStart, is->mEnd, is->mScore));
        }
    }
    // sort
    for(int32_t i = 0; i < tid; ++i)  std::sort(oci[i].begin(), oci[i].end());
}

void VCFMerger::writeIntervals(ContigIntervals& ci, ContigMap& cmap, int32_t svtin){
    // open output bcf file
    htsFile* fp = hts_open(outfile.c_str(), "wb");
    bcf_hdr_t* hdrOut = bcf_hdr_init("w");
    // write vcf header
    bcf_hdr_append(hdrOut, "##ALT=<ID=DEL,Description=\"Deletion\">");
    bcf_hdr_append(hdrOut, "##ALT=<ID=DUP,Description=\"Duplication\">");
    bcf_hdr_append(hdrOut, "##ALT=<ID=INV,Description=\"Inversion\">");
    bcf_hdr_append(hdrOut, "##ALT=<ID=BND,Description=\"Translocation\">");
    bcf_hdr_append(hdrOut, "##ALT=<ID=INS,Description=\"Insertion\">");
    bcf_hdr_append(hdrOut, "##FILTER=<ID=LowQual,Description=\"Poor quality and insufficient number of PEs and SRs.\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=PEMAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of split-reads\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=SRALNQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Predicted length of the insertion\">");
    bcf_hdr_append(hdrOut, "##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Predicted microhomology length using a max. edit distance of 2\">");
    // append reference contig information
    std::map<int32_t, std::string> imap;
    for(auto& e: cmap) imap[e.second] = e.first;
    for(auto& e: imap){
        std::string ctgstr = "##contig=<ID=" + e.second + ">";
        bcf_hdr_append(hdrOut, ctgstr.c_str());
    }
    bcf_hdr_add_sample(hdrOut, NULL);
    bcf_hdr_write(fp, hdrOut);
    // duplicate filter
    int32_t numseq = imap.size();
    std::vector<std::set<std::pair<int32_t, int32_t>>> gis(numseq);
    // parse input bcf files
    bcf1_t* rout = bcf_init();
    std::vector<htsFile*> ifp(infiles.size());
    std::vector<bcf_hdr_t*> hdr(infiles.size());
    std::vector<bcf1_t*> rec(infiles.size());
    std::vector<bool> eof(infiles.size());
    uint32_t allEOF = 0;
    for(uint32_t c = 0; c < infiles.size(); ++c){
        ifp[c] = hts_open(infiles[c].c_str(), "r");
        hdr[c] = bcf_hdr_read(ifp[c]);
        rec[c] = bcf_init();
        if(bcf_read(ifp[c], hdr[c], rec[c]) == 0){
            bcf_unpack(rec[c], BCF_UN_INFO);
            eof[c] = false;
        }else{
            ++allEOF;
            eof[c] = true;
        }
    }
    int32_t nsvend = 0;
    int32_t* svend = NULL;
    int32_t npe = 0;
    int32_t* pe = NULL;
    int32_t nsr = 0;
    int32_t* sr = NULL;
    int32_t ninslen = 0;
    int32_t* inslen = NULL;
    int32_t nhomlen = 0;
    int32_t* homlen = NULL;
    int32_t npemapq = 0;
    int32_t* pemapq = NULL;
    int32_t nsrmapq = 0;
    int32_t* srmapq = NULL;
    int32_t nsralnq = 0;
    float* sralnq = NULL;
    int32_t nct = 0;
    char* ct = NULL;
    int32_t nsvt = 0;
    char* svt = NULL;
    int32_t nchr2 = 0;
    char* chr2 = NULL;
    int32_t ncipos = 0;
    int32_t* cipos = NULL;
    int32_t nciend = 0;
    int32_t* ciend = NULL;
    int32_t ncons = 0;
    char* cons = NULL;
    while(allEOF < infiles.size()){
        // find next sorted record(by pos and tid)
        int32_t idx = -1;
        for(uint32_t c = 0; c < infiles.size(); ++c){
            if(!eof[c]){
                if((idx < 0) || rec[idx]->rid > rec[c]->rid || (rec[idx]->rid == rec[c]->rid && rec[idx]->pos > rec[c]->pos)){
                    idx = c;
                }
            }
        }
        // get SV type
        int32_t recsvt = -1;
        if(bcf_get_info_string(hdr[idx], rec[idx], "SVTYPE", &svt, &nsvt) > 0 &&
           bcf_get_info_string(hdr[idx], rec[idx], "CT", &ct, &nct) > 0){
            recsvt = svutil::str2svt(ct, svt);
        }
        if(recsvt == svtin){
            // check PASS
            bool passFilter = true;
            if(filterForPass) passFilter = (bcf_has_filter(hdr[idx], rec[idx], const_cast<char*>("PASS")) == 1);
            // check PRECISE
            bool precise = false, passPrecise = true;
            if(bcf_get_info_flag(hdr[idx], rec[idx], "PRECISE", 0, 0) > 0) precise = true;
            if(filterForPrecise && (!precise)) passPrecise = false;
            // check PASS and PRECISE
            if(passPrecise && passFilter){
                // correct sv size
                std::string chr1Name = bcf_hdr_id2name(hdr[idx], rec[idx]->rid);
                int32_t tid = cmap[chr1Name];
                int32_t svStart = rec[idx]->pos;
                int32_t svEnd = svStart + 1;
                if(bcf_get_info_int32(hdr[idx], rec[idx], "END", &svend, &nsvend) > 0) svEnd = *svend;
                int32_t inslenVal = 0;
                if(bcf_get_info_int32(hdr[idx], rec[idx], "INSLEN", &inslen, &ninslen) > 0) inslenVal = *inslen;
                // Parse INFO fields
                if(std::strcmp(svt, "BND") == 0 || (std::strcmp(svt, "INS") == 0 && inslenVal >= minsize) ||
                   (std::strcmp(svt, "INS") != 0 && std::strcmp(svt, "BND") != 0 && svEnd - svStart > minsize && svEnd - svStart < maxsize)){
                    int32_t peSupport = 0, srSupport = 0; // PE SR
                    int32_t peMapQual = 0, srMapQual = 0; // PEMAPQ SRMAPQ
                    if(bcf_get_info_int32(hdr[idx], rec[idx], "PE", &pe, &npe) > 0) peSupport = *pe;
                    if(bcf_get_info_int32(hdr[idx], rec[idx], "SR", &sr, &nsr) > 0) srSupport = *sr;
                    if(bcf_get_info_int32(hdr[idx], rec[idx], "PEMAPQ", &pemapq, &npemapq) > 0) peMapQual = *pemapq;
                    if(bcf_get_info_int32(hdr[idx], rec[idx], "SRMAPQ", &srmapq, &nsrmapq) > 0) srMapQual = *srmapq;
                    std::string chr2Name = chr1Name; // CHR2
                    if(bcf_get_info_string(hdr[idx], rec[idx], "CHR2", &chr2, &nchr2) > 0) chr2Name = chr2;
                    // get quality for this SV
                    int32_t svScore = 3 * srSupport * srMapQual + peSupport * peMapQual;
                    if(bcf_hdr_id2int(hdr[idx], BCF_DT_ID, "SCORE") >= 0){
                        int32_t nscore = 0;
                        if(bcf_hdr_id2type(hdr[idx], BCF_HL_INFO, bcf_hdr_id2int(hdr[idx], BCF_DT_ID, "SCORE")) == BCF_HT_INT){
                            int32_t* score = NULL;
                            bcf_get_info_int32(hdr[idx], rec[idx], "SCORE", &score, &nscore);
                            svScore = *score;
                            free(score);
                        }else if(bcf_hdr_id2type(hdr[idx], BCF_HL_INFO, bcf_hdr_id2int(hdr[idx], BCF_DT_ID, "SCORE")) == BCF_HT_REAL){
                            float* score = NULL;
                            bcf_get_info_float(hdr[idx], rec[idx], "SCORE", & score, &nscore);
                            svScore = 10000 * *score;
                            free(score);
                        }
                    }
                    // check this interval is selected
                    Interval ti(svStart, svEnd, svScore);
                    auto iter = std::lower_bound(ci[tid].begin(), ci[tid].end(), Interval(svStart, svEnd, svScore));
                    bool validInterval = false;
                    for(; iter != ci[tid].end() && iter->mStart == svStart; ++iter){
                        if(iter->mStart == svStart && iter->mEnd == svEnd && iter->mScore == svScore){
                            if(gis[tid].find({svStart, svEnd}) == gis[tid].end()){
                                validInterval = true;
                                gis[tid].insert({svStart, svEnd});
                            }
                            break;
                        }
                    }
                    if(validInterval){
                        // Fetch missing INFO fields
                        int32_t homlenVal = 0;
                        if(bcf_get_info_int32(hdr[idx], rec[idx], "HOMLEN", &homlen, &nhomlen) > 0) homlenVal = *homlen;
                        bcf_get_info_int32(hdr[idx], rec[idx], "CIPOS", &cipos, &ncipos);
                        bcf_get_info_int32(hdr[idx], rec[idx], "CIEND", &ciend, &nciend);
                        float srAlnQual = 0;
                        if(bcf_get_info_float(hdr[idx], rec[idx], "SRALNQ", &sralnq, &nsralnq) > 0) srAlnQual = *sralnq;
                        std::string cs;
                        if(precise){
                            bcf_get_info_string(hdr[idx], rec[idx], "CONSENSUS", &cons, & ncons);
                            cs = cons;
                        }
                        // create new record
                        rout->rid = bcf_hdr_name2id(hdrOut, chr1Name.c_str());
                        rout->pos = rec[idx]->pos;
                        rout->qual = svScore;
                        std::string id = svutil::addID(svtin);
                        id.append(std::to_string(svcounter++));
                        bcf_update_id(hdrOut, rout, id.c_str());
                        std::string ref = rec[idx]->d.allele[0];
                        std::string alt = rec[idx]->d.allele[1];
                        std::string alleles = ref + "," + alt;
                        bcf_update_alleles_str(hdrOut, rout, alleles.c_str());
                        int32_t tmpPass = bcf_hdr_id2int(hdrOut, BCF_DT_ID, "PASS");
                        bcf_update_filter(hdrOut, rout, &tmpPass, 1);
                        // add INFO fields
                        if(precise) bcf_update_info_flag(hdrOut, rout, "PRECISE", NULL, 1);
                        else bcf_update_info_flag(hdrOut, rout, "IMPRECISE", NULL, 1);
                        bcf_update_info_string(hdrOut, rout, "SVTYPE", svutil::addID(svtin).c_str());
                        bcf_update_info_string(hdrOut, rout, "CHR2", chr2Name.c_str());
                        bcf_update_info_int32(hdrOut, rout, "END", &svEnd, 1);
                        bcf_update_info_int32(hdrOut, rout, "PE", &peSupport, 1);
                        int32_t tmpi = peMapQual;
                        bcf_update_info_int32(hdrOut, rout, "EAMAPQ", &tmpi, 1);
                        bcf_update_info_string(hdrOut, rout, "CT", svutil::addOrientation(svtin).c_str());
                        bcf_update_info_int32(hdrOut, rout, "CIPOS", cipos, 2);
                        bcf_update_info_int32(hdrOut, rout, "CIEND", ciend, 2);
                        if(precise){
                            bcf_update_info_int32(hdrOut, rout, "INSLEN", &inslenVal, 1);
                            bcf_update_info_int32(hdrOut, rout, "HOMLEN", &homlenVal, 1);
                            bcf_update_info_int32(hdrOut, rout, "SR", &srSupport, 1);
                            tmpi = srMapQual;
                            bcf_update_info_int32(hdrOut, rout, "SRMAPQ", &tmpi, 1);
                            bcf_update_info_float(hdrOut, rout, "SRALNQ", &srAlnQual, 1);
                            bcf_update_info_string(hdrOut, rout, "CONSENSUS", cs.c_str());
                        }
                        // write record
                        bcf_write(fp, hdrOut, rout);
                        bcf_clear1(rout);
                    }
                }
            }
        }
        if(bcf_read(ifp[idx], hdr[idx], rec[idx]) == 0) bcf_unpack(rec[idx], BCF_UN_INFO);
        else{
            ++allEOF;
            eof[idx] = true;
        }
    }
    // free memory
    if(svend) free(svend);
    if(pe) free(pe);
    if(sr) free(sr);
    if(inslen) free(inslen);
    if(homlen) free(homlen);
    if(pemapq) free(pemapq);
    if(srmapq) free(srmapq);
    if(sralnq) free(sralnq);
    if(ct) free(ct);
    if(svt) free(svt);
    if(chr2) free(chr2);
    if(cipos) free(cipos);
    if(ciend) free(ciend);
    if(cons) free(cons);
    // clean-up
    for(uint32_t c = 0; c < infiles.size(); ++c){
        bcf_hdr_destroy(hdr[c]);
        bcf_destroy(rec[c]);
        bcf_close(ifp[c]);
    }
    // close output bcf file
    bcf_destroy(rout);
    bcf_hdr_destroy(hdrOut);
    hts_close(fp);
    // Build index
    bcf_index_build(outfile.c_str(), 14);
}

void VCFMerger::mergeBCFs(std::vector<std::string>& ifiles){
    // parse input bcf files
    std::vector<htsFile*> ifp(ifiles.size());
    std::vector<bcf_hdr_t*> hdr(ifiles.size());
    std::vector<bcf1_t*> rec(ifiles.size());
    std::vector<bool> eof(ifiles.size());
    uint32_t allEOF = 0;
    for(uint32_t c = 0; c < ifiles.size(); ++c){
        ifp[c] = bcf_open(ifiles[c].c_str(), "r");
        hdr[c] = bcf_hdr_read(ifp[c]);
        rec[c] = bcf_init();
        if(bcf_read(ifp[c], hdr[c], rec[c]) == 0) eof[c] = false;
        else{
            eof[c] = true;
            ++allEOF;
        }
    }
    // open output bcf file
    htsFile* fp = bcf_open(outfile.c_str(), "wb");
    bcf_hdr_t* hdrOut = bcf_hdr_dup(hdr[0]);
    bcf_hdr_write(fp, hdrOut);
    // merge files
    while(allEOF < ifiles.size()){
        // find next sorted record(by tid and pos)
        int32_t idx = -1;
        for(uint32_t c = 0; c < ifiles.size(); ++c){
            if(!eof[c]){
                if((idx < 0) || (rec[idx]->rid > rec[c]->rid) || (rec[idx]->rid == rec[c]->rid && rec[idx]->pos > rec[c]->pos)) idx = c;
            }
        }
        // write record
        bcf_write(fp, hdrOut, rec[idx]);
        // fetch next record
        if(bcf_read(ifp[idx], hdr[idx], rec[idx])){
            ++allEOF;
            eof[idx] = true;
        }
    }
    // clean-up
    for(uint32_t c = 0; c < ifiles.size(); ++c){
        bcf_hdr_destroy(hdr[c]);
        bcf_destroy(rec[c]);
        bcf_close(ifp[c]);
    }
    // close bcf file
    bcf_hdr_destroy(hdrOut);
    bcf_close(fp);
    // build index
    bcf_index_build(outfile.c_str(), 14);
}

void VCFMerger::mergeOneType(int32_t svt){
    util::loginfo("Beg Processing SVT: " + std::to_string(svt));
    // get contig ~ tid map
    ContigMap cmap;
    getContigMap(cmap);
    // get contig intervals
    ContigIntervals ici(cmap.size());
    getIntervals(ici, cmap, svt);
    int32_t svc = 0;
    for(auto& e: ici) svc += e.size();
    util::loginfo("Got " + std::to_string(svc) + " raw Intervals for SVT: " + std::to_string(svt));
    // merge intervals
    ContigIntervals oci(cmap.size());
    mergeIntervals(ici, oci, svt);
    svc = 0;
    for(auto& e: oci) svc += e.size();
    util::loginfo("Got " + std::to_string(svc) + " merged Intervals for SVT: " + std::to_string(svt));
    // write SVs in merged intervals
    writeIntervals(oci, cmap, svt);
    util::loginfo("End processing SVT: " + std::to_string(svt));
}


void VCFMerger::mergeAllType(){
    // merge all files
    util::loginfo("Beg merging each type of SVs");
    std::string origOutfile = outfile;
    std::vector<std::string> svtOutPath(9);
    for(int32_t svt = 0; svt < 9; ++svt){
        svtOutPath[svt] = util::joinpath(tmpdir, "tmp." + std::to_string(svt) + ".bcf");
        if(infiles.size() < chunksize){// merge in one go
            outfile = svtOutPath[svt];
            mergeOneType(svt);
        }else{// merge in chunks
            std::vector<std::string> origInFiles = infiles;
            int32_t chunks = (infiles.size() - 1) / chunksize + 1;
            std::vector<std::string> chunkOutPath(chunks);
            for(int32_t ic = 0; ic < chunks; ++ic){
                chunkOutPath[ic] = util::joinpath(tmpdir, "chk.tmp." + std::to_string(ic) + ".bcf");
                infiles.clear();
                for(uint32_t k = ic * chunksize; k < (ic + 1) * chunksize && k < origInFiles.size(); ++k){
                    infiles.push_back(origInFiles[k]);
                }
                outfile = chunkOutPath[ic];
                mergeOneType(svt);
            }
            // merge chunks
            infiles = chunkOutPath;
            outfile = svtOutPath[svt];
            // reset vaf and coverage because these are site lists!
            float origVaf = vaf;
            int32_t origCov = coverage;
            vaf = 0;
            coverage = 0;
            mergeOneType(svt);
            vaf = origVaf;
            coverage = origCov;
            // clean-up
            for(int32_t ic = 0; ic < chunks; ++ic){
                remove(chunkOutPath[ic].c_str());
                std::string bcfidx = chunkOutPath[ic] + ".csi";
                remove(bcfidx.c_str());
            }
            infiles = origInFiles;
        }
    }
    util::loginfo("End merging each type of SVs");
    // Merge temporary files
    outfile = origOutfile;
    util::loginfo("Beg merging all types of SVs");
    mergeBCFs(svtOutPath);
    util::loginfo("End merging all types of SVs");
    for(int32_t svt = 0; svt < 9; ++svt){
        remove(svtOutPath[svt].c_str());
        std::string bcfidx = svtOutPath[svt] + ".csi";
        remove(bcfidx.c_str());
    }
    // Remove tmpdir
    rmdir(tmpdir.c_str());
}

int main(int argc, char** argv){
    if(argc == 1){
        std::string helpCMD = std::string(argv[0]) + " -h";
        std::system(helpCMD.c_str());
        return 0;
    }
    // parse commandline arguments
    VCFMerger* m = new VCFMerger();
    CLI::App app("program: " + std::string(argv[0]) + "\n" + m->softEnv->cmp);
    app.get_formatter()->column_width(36);
    app.add_option("-i,--in", m->infilelist, "input bcf file list to be merged")->required(true)->check(CLI::ExistingFile);
    app.add_option("-o,--out", m->outfile, "merged bcf file output path", true);
    app.add_option("-t,--tmp", m->tmpdir, "temp directory to store intermediate files", true);
    app.add_option("-v,--vaf", m->vaf, "min VAF for SV to be merged", true);
    app.add_option("-c,--cov", m->coverage, "min depth for SV to be merged", true);
    app.add_option("--covratio", m->recoverlap, "min overlap ratio needed for two dup SVs", true);
    app.add_option("--bpoffset", m->bpoffset, "max position offset allowed for two dup SV", true);
    app.add_option("--minsize", m->minsize, "min size of SV to be merged", true);
    app.add_option("--maxsize", m->maxsize, "max size of SV to be merged", true);
    app.add_option("--chunksize", m->chunksize, "max input files to be merged at one batch", true);
    app.add_flag("--passfilter", m->filterForPass, "keep PASS records only if set");
    app.add_flag("--preciseFilter", m->filterForPrecise, "keep PRECISE records only if set");
    CLI_PARSE(app, argc, argv);
    m->update(argc, argv);
    m->mergeAllType();
    delete m;
}
