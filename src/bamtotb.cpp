#include "bamtotb.h"

void BamToTable::b2r(bam1_t* b, bam_hdr_t* h, BamRec& br, int32_t id){
    br.chr = h->target_name[b->core.tid];
    br.pos = b->core.pos;
    br.mchr = h->target_name[b->core.mtid];
    br.mpos = b->core.mpos;
    br.cigar = bamutil::getCigar(b);
    br.mcigar = bamutil::getStrTag(b, "MC");
    br.sa = bamutil::getStrTag(b, "SA");
    br.barcode = bamutil::getStrTag(b, "BC");
    br.seq = bamutil::getSeq(b);
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
                br.rbp = refpos;
                if(i == 0) lsc = true;
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
        for(auto& e: pcigar){
            char opchr = e.second;
            int oplen = e.first;
            if(opchr == 'M' || opchr == '=' || opchr == 'X' || opchr == 'D' || opchr == BAM_CREF_SKIP){
                refpos += oplen;
            }else if(opchr == 'S'){
                br.sbp = refpos;
                break;
            }
        }
    }else{
        br.sbp = -1;
    }
    br.svid = id;
    br.qname = bamutil::getQName(b);
    br.svrt = bamutil::getIntTag(b, "ST");
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
    // fetch svids
    std::map<int32_t, FRExtraInfo> fim;
    getFRExtraInfo(fim);
    // fetch bam records
    BamRecVector brecs;
    samFile* fp = sam_open(svbam.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    while(sam_read1(fp, h, b) >= 0){
        uint8_t* data = bam_aux_get(b, "ZF");
        if(data){
            int32_t id= bam_aux2i(data);
            auto iter = fim.find(id);
            if(iter != fim.end()){
                BamRec br;
                b2r(b, h, br, id);
                br.fsgene = iter->second.fsgene;
                brecs.push_back(br);
            }
        }
    }
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
        lxw_worksheet* rsheet = workbook_add_worksheet(workbook, "FusionReads");
        int ttl = lxwutil::lines2sheet(rsheet, oss.str(), format);
        worksheet_autofilter(rsheet, 0, 0, std::max(0, ttl - 2), 0);
        workbook_close(workbook);
    }
    // release memory
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
}
