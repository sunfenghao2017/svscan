#include "bamtotb.h"

void BamToTable::lines2sheet(lxw_worksheet* sheet, const std::string& buf, lxw_format* fmt){
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
}

std::string BamToTable::b2s(bam1_t* b, bam_hdr_t* h){
    std::stringstream ss;
    ss << h->target_name[b->core.tid] << "\t" << b->core.pos << "\t"; // chr pos 
    ss << h->target_name[b->core.mtid] << "\t" << b->core.mpos << "\t"; // mchr mpos
    ss << bamutil::getCigar(b) << "\t"; // cigar
    std::string mc = bamutil::getStrTag(b, "MC"); // mcigar
    if(mc.empty()) ss << "-\t";
    else ss << mc << "\t";
    std::string sa = bamutil::getStrTag(b, "SA"); // sa
    if(sa.empty()) ss << "-\t";
    else ss << sa << "\t";
    ss << bamutil::getSeq(b) << "\n"; // seq
    return ss.str();
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
    if(!sstsv.empty()){
        std::ifstream fr(sstsv);
        std::string line;
        std::getline(fr, line);
        std::vector<std::string> vstr;
        while(std::getline(fr, line)){
            util::split(line, vstr, "\t");
            svids.insert(std::atoi(vstr[svidf].c_str()));
        }
        fr.close();
    }
}

void BamToTable::b2t(){
    std::set<int32_t> svids;
    getsvid(svids);
    lxw_workbook* workbook = new_workbook(bamtb.c_str());
    lxw_format* format = workbook_add_format(workbook);
    format_set_align(format, LXW_ALIGN_LEFT);
    format_set_align(format, LXW_ALIGN_VERTICAL_BOTTOM);
    std::string header = "chr\tpos\tmchr\tmpos\tcigar\tmcigar\tsa\tseq\n";
    std::map<int32_t, std::pair<std::string, lxw_worksheet*>> svid2result;
    for(auto& id: svids){
        svid2result[id].first.append(header);
        svid2result[id].second = workbook_add_worksheet(workbook, std::to_string(id).c_str());
    }
    samFile* fp = sam_open(svbam.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    while(sam_read1(fp, h, b) >= 0){
        uint8_t* data = bam_aux_get(b, "ZF");
        if(data){
            int32_t id= bam_aux2i(data);
            auto iter = svid2result.find(id);
            if(iter != svid2result.end()){
                iter->second.first.append(b2s(b, h));
            }
        }
    }
    for(auto iter = svid2result.begin(); iter != svid2result.end(); ++iter){
        lines2sheet(iter->second.second, iter->second.first, format);
    }
    workbook_close(workbook);
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
}
