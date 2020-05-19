#include "bamtotb.h"

void BamToTable::b2r(bam1_t* b, bam_hdr_t* h, BamRec& br, int32_t id){
    br.chr = h->target_name[b->core.tid];
    br.pos = b->core.pos + 1;
    br.cigar = bamutil::getCigar(b);
    br.svrt = bamutil::getIntTag(b, "ST");
    if(b->core.mtid >= 0 && (!(b->core.flag & BAM_FMUNMAP))){
        br.mchr = h->target_name[b->core.mtid];
        br.mpos = b->core.mpos + 1;
        br.mcigar = bamutil::getStrTag(b, "MC");
    }
    br.sa = bamutil::getStrTag(b, "SA");
    br.barcode = bamutil::getStrTag(b, "BC");
    br.seq = bamutil::getSeq(b);
    if(br.svrt == 1){
        br.lseq = br.seq;
        br.rbp = -1;
    }else{
        std::pair<int32_t, int32_t> scl = bamutil::getSoftClipLength(b);
        if(scl.first + scl.second == 0){
            br.lseq = br.seq;
            br.rbp = -1;
        }else if((scl.first > 0) ^ (scl.second > 0)){
            if(scl.first){
                br.lseq = br.seq.substr(0, scl.first);
                br.tseq = br.seq.substr(scl.first);
            }
            if(scl.second){
                br.lseq = br.seq.substr(0, br.seq.length() - scl.second);
                br.tseq = br.seq.substr(br.seq.length() - scl.second);
            }
            uint32_t* cigar = bam_get_cigar(b);
            int refpos = b->core.pos;
            bool lsc = false;
            for(uint32_t i = 0; i < b->core.n_cigar; ++i){
                int opint = bam_cigar_op(cigar[i]);
                int oplen = bam_cigar_oplen(cigar[i]);
                if(opint == BAM_CMATCH || opint == BAM_CEQUAL || opint == BAM_CDIFF || opint == BAM_CDEL || opint == BAM_CREF_SKIP){
                    refpos += oplen;
                }else if(opint == BAM_CSOFT_CLIP){
                    if(i == 0){
                        br.rbp = refpos + 1;
                        lsc = true;
                    }else{
                        br.rbp = refpos;
                    }
                    break;
                }
            }
            if(lsc){
                br.lhit = bamutil::getIntTag(b, "SH");
                br.thit = bamutil::getIntTag(b, "PH");
            }else{
                br.lhit = bamutil::getIntTag(b, "PH");
                br.thit = bamutil::getIntTag(b, "SH");
            }
        }
    }
    if(br.sa.size()){
        // get optimal SA
        std::vector<std::string> cvs;
        std::vector<std::string> vstr;
        util::split(br.sa, cvs, ";");
        if(cvs[1].empty()){
            br.sa = cvs[0];
        }else{
            for(uint32_t cvidx = 0; cvidx < cvs.size() - 1; ++cvidx){
                util::split(cvs[cvidx], vstr, ",");
                if(vstr[3].find_first_of("SH") == vstr[3].find_last_of("SH")){
                    br.sa = cvs[cvidx];
                    break;
                }
            }
        }
        util::split(br.sa, vstr, ",");
        int32_t refpos = std::atoi(vstr[1].c_str());
        std::vector<std::pair<int32_t, char>> pcigar;
        bamutil::parseCigar(vstr[3], pcigar);
        for(uint32_t si = 0; si < pcigar.size(); ++si){
            char opchr = pcigar[si].second;
            int oplen = pcigar[si].first;
            if(opchr == 'M' || opchr == '=' || opchr == 'X' || opchr == 'D' || opchr == BAM_CREF_SKIP){
                refpos += oplen;
            }else if(opchr == 'S'){
                if(si == 0){
                    br.sbp = refpos;
                }else{
                    br.sbp = refpos-1;
                }
                break;
            }
        }
    }else{
        br.sbp = -1;
    }
    br.svid = id;
    br.qname = bamutil::getQName(b);
    br.read1 = (b->core.flag & BAM_FREAD1);
    if(b->core.flag & BAM_FREVERSE) br.strand = '-';
    else br.strand = '+';
    if(b->core.flag & BAM_FMREVERSE) br.mstrand = '-';
    else br.mstrand = '+';
}

void BamToTable::getFRExtraInfo(std::map<int32_t, FRExtraInfo>& fim){
    if(!fstsv.empty()){
        std::ifstream fr(fstsv);
        std::string line;
        std::getline(fr, line);
        std::vector<std::string> vstr;
        while(std::getline(fr, line)){
            util::split(line, vstr, "\t");
            FRExtraInfo fi;
            fi.svid = std::atoi(vstr[svidf].c_str());
            fi.fsgene = vstr[fsidf];
            fim[fi.svid] = fi;
        }
        fr.close();
    }

    for(uint32_t i = 0; i < usrid.size(); ++i){
        FRExtraInfo fi;
        fi.svid = usrid[i];
        if(usrid.size() == fsgene.size()){
            fi.fsgene = fsgene[i];
        }else{
            fi.fsgene = "-";
        }
        fim[usrid[i]] = fi;
    }
}

void BamToTable::b2t(){
    if(bamtt.empty() && bamtb.empty()) return; // return if both types of output does not needed
    regs = new BedRegs();
    regs->mCR = cr_init();
    // fetch svids
    std::map<int32_t, FRExtraInfo> fim;
    getFRExtraInfo(fim);
    // fetch all bam records
    BamRecVector brecs;
    samFile* fp = sam_open(svbam.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    std::vector<bam1_t*> bamrecs;
    bam1_t* b = bam_init1();
    while(sam_read1(fp, h, b) >= 0){
        uint8_t* data = bam_aux_get(b, "ZF");
        if(data){
            int32_t id = bam_aux2i(data);
            auto iter = fim.find(id);
            if(iter != fim.end()){
                BamRec br;
                b2r(b, h, br, id);
                br.fsgene = iter->second.fsgene;
                brecs.push_back(br);
                uint8_t *svtd = bam_aux_get(b, "ST");
                if(svtd){
                    int svtv = bam_aux2i(svtd);
                    if(svtv){
                        int64_t mask = id;
                        mask = mask << 32;
                        if(b->core.flag & BAM_FREAD1) mask |= BAM_FREAD2;
                        else mask |= BAM_FREAD1;
                        r2get[bam_get_qname(b)] = mask;
                        if(!(b->core.flag & BAM_FUNMAP)){
                            cr_add(regs->mCR, sam_hdr_tid2name(h, b->core.mtid), b->core.mpos, b->core.mpos + 2, 0);
                        }
                    }
                }
                bamrecs.push_back(b);
                b = bam_init1();
            }
        }
    }
    sam_close(fp);
    // fetch other pe supporting reads
    fp = sam_open(ttbam.c_str(), "r");
    hts_idx_t *idx = sam_index_load(fp, ttbam.c_str());
    cr_index2(regs->mCR, 1);
    cgranges_t *cr = regs->mCR;
    for(int32_t ctg_id = 0; ctg_id < cr->n_ctg; ++ctg_id){
        int64_t i, *xb = 0, max_b = 0, n = 0;
        n = cr_overlap_int(cr, ctg_id, 0, INT_MAX, &xb, &max_b);
        for(i = 0; i < n; ++i){
            const char* chr = cr->ctg[ctg_id].name;
            int32_t beg = cr_start(cr, xb[i]), end = cr_end(cr, xb[i]);
            hts_itr_t *itr = sam_itr_queryi(idx, sam_hdr_name2tid(h, chr), beg, end);
            while(sam_itr_next(fp, itr, b) >= 0){
                auto iter = r2get.find(bam_get_qname(b));
                if(iter != r2get.end()){
                    if(b->core.flag & (iter->second & 0xffff)){
                        int32_t id = iter->second >> 32;
                        bam_aux_update_int(b, "SZ", id);
                        bam_aux_update_int(b, "ST", 1);
                        for(uint32_t k = 0; k < brecs.size(); ++k){
                            if(brecs[k].svid == id && brecs[k].svrt == 1 && brecs[k].qname == bam_get_qname(b)){
                                brecs[k].tseq = bamutil::getSeq(b);
                                brecs[k].mcigar = bamutil::getCigar(b);
                                bamrecs.push_back(b);
                                b = bam_init1();
                                break;
                            }
                        }
                    }
                }
            }
        }
        free(xb);
    }
    sam_close(fp);
    bam_destroy1(b);
    // sort bam
    samFile* of = sam_open(newbam.c_str(), "wb");
    assert(sam_hdr_write(of, h) >= 0);
    std::sort(bamrecs.begin(), bamrecs.end(), BamComp());
    for(auto& e: bamrecs){
        assert(sam_write1(of, h, e) >= 0);
        bam_destroy1(e);
    }
    sam_close(of);
    assert(sam_index_build(newbam.c_str(), 0) == 0);
    std::sort(brecs.begin(), brecs.end());
    // bam records to string buffer
    std::stringstream oss;
    oss << BamRec::getHeader();
    for(uint32_t i = 0; i < brecs.size(); ++i){
        oss << brecs[i].toStr();
    }
    // store into txt if needed
    if(!bamtt.empty()){
        std::ofstream fw(bamtt);
        fw << oss.str();
        fw.close();
    }
    // store into excel if needed
    if(!bamtb.empty()){
        lxw_workbook* workbook = workbook_new(bamtb.c_str());
        lxw_format* format = workbook_add_format(workbook);
        format_set_align(format, LXW_ALIGN_LEFT);
        format_set_align(format, LXW_ALIGN_VERTICAL_BOTTOM);
        lxw_worksheet* rsheet = workbook_add_worksheet(workbook, "FusionMols");
        int ttl = lxwutil::lines2sheet(rsheet, oss.str(), format);
        worksheet_autofilter(rsheet, 0, 0, std::max(0, ttl - 2), 0);
        workbook_close(workbook);
    }
}
