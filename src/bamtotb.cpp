#include "bamtotb.h"

int BamToTable::lines2sheet(lxw_worksheet* sheet, const std::string& buf, lxw_format* fmt){
    int row = 0, col = 0;
    std::vector<size_t> vclen;
    std::vector<std::string> vline;
    util::split(buf, vline, "\n");
    std::vector<std::string> vrec;
    size_t maxcol = 0;
    for(auto& line: vline){
        util::split(line, vrec, "\t");
        maxcol = std::max(maxcol, vrec.size());
    }
    vclen.resize(maxcol, 0);
    for(auto& line: vline){
         if(line.empty()){
             ++row;
             continue;
         }
         util::split(line, vrec, "\t");
         col = 0;
         for(size_t i = 0; i < vrec.size(); ++i){
            vclen[i] = std::max(vclen[i], vrec[i].length());
            char* p = NULL;
            float val = std::strtof(vrec[i].c_str(), &p);
            if(*p) worksheet_write_string(sheet, row, col++, vrec[i].c_str(), NULL);
            else worksheet_write_number(sheet, row, col++, val, NULL);
         }
         ++row;
     }
     for(size_t i = 0; i < vclen.size(); ++i) worksheet_set_column(sheet, i, i, vclen[i], fmt);
     return row;
}

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
    if(scl.first + scl.second == 0) br.lseq = br.seq;
    else if((scl.first > 0) ^ (scl.second > 0)){
        if(scl.first){
            br.lseq = br.seq.substr(0, scl.first);
            br.tseq = br.seq.substr(scl.first);
        }
        if(scl.second){
            br.lseq = br.seq.substr(0, br.seq.length() - scl.second);
            br.tseq = br.seq.substr(br.seq.length() - scl.second);
        }
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

void BamToTable::getsvid(std::set<int32_t>& svids){
    if(!fstsv.empty()){
        std::ifstream fr(fstsv);
        std::string line;
        std::getline(fr, line);
        std::vector<std::string> vstr;
        while(std::getline(fr, line)){
            util::split(line, vstr, "\t");
            svids.insert(std::atoi(vstr[svidf].c_str()));
        }
        fr.close();
    }
    for(auto& e: usrid){
        svids.insert(e);
    }
}

void BamToTable::b2t(){
    std::set<int32_t> svids;
    getsvid(svids);
    lxw_workbook* workbook = workbook_new(bamtb.c_str());
    lxw_format* format = workbook_add_format(workbook);
    format_set_align(format, LXW_ALIGN_LEFT);
    format_set_align(format, LXW_ALIGN_VERTICAL_BOTTOM);
    BamRecVector brecs;
    samFile* fp = sam_open(svbam.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    while(sam_read1(fp, h, b) >= 0){
        uint8_t* data = bam_aux_get(b, "ZF");
        if(data){
            int32_t id= bam_aux2i(data);
            auto iter = svids.find(id);
            if(iter != svids.end()){
                BamRec br;
                b2r(b, h, br, id);
                brecs.push_back(br);
            }
        }
    }
    std::sort(brecs.begin(), brecs.end());
    std::stringstream oss;
    oss << BamRec::getHeader();
    for(uint32_t i = 0; i < brecs.size(); ++i){
        oss << brecs[i].toStr();
    }
    lxw_worksheet* rsheet = workbook_add_worksheet(workbook, "FusionReads");
    int ttl = lines2sheet(rsheet, oss.str(), format);
    worksheet_autofilter(rsheet, 0, 0, std::max(0, ttl - 2), 0);
    workbook_close(workbook);
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
}
