#include "svdebug.h"

void SVDebug::getReg(const std::string& trs, tbx_t* tbx, htsFile* tfp, GeneRegion& reg){
    hts_itr_t* itr = tbx_itr_querys(tbx, trs.c_str());
    kstring_t rec = {0, 0, 0};
    std::vector<std::string> vstr;
    std::vector<GeneRegion> regs;
    while(tbx_itr_next(tfp, tbx, itr, &rec) >= 0){
        util::split(rec.s, vstr, "\t");
        if(vstr[3] == "exon"){
            GeneRegion gr;
            gr.beg = vstr[7];
            gr.end = vstr[8];
            gr.chr = vstr[6];
            gr.strand = vstr[9];
            regs.push_back(gr);
        }
    }
    hts_itr_destroy(itr);
    if(regs.size() > 0){
        reg.strand = regs[0].strand;
        reg.chr = regs[0].chr;
        if(reg.strand[0] == '-'){
            reg.end = regs[0].end;
            reg.beg = regs[regs.size() - 1].beg;
        }else{
            reg.end = regs[regs.size() - 1].end;
            reg.beg = regs[0].beg;
        }
    }else{
        util::errorExit("Transcript : " + trs + "has no exons");
    }
}

void SVDebug::debugDNA(){
    if(annodb.empty()) util::errorExit("annodb must be provided in DNA SV debugging");
    // fetch gene coordinates
    htsFile* tfp = hts_open(annodb.c_str(), "r");
    tbx_t* tbx = tbx_index_load(annodb.c_str());
    std::map<std::string, std::string> g2tmap;
    util::makeMapPairFromFileByLine(gene2trans, g2tmap);
    std::string htrs = g2tmap[hgene];
    std::string ttrs = g2tmap[tgene];
    GeneRegion hreg, treg;
    getReg(htrs, tbx, tfp, hreg);
    getReg(ttrs, tbx, tfp, treg);
    std::string hqreg = hreg.toString();
    std::string tqreg = treg.toString();
    util::loginfo(hgene + ": " + hqreg);
    util::loginfo(tgene + ": " + tqreg);
    hts_close(tfp);
    tbx_destroy(tbx);
    // open bam and idx
    samFile* sfp = sam_open(inbam.c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(sfp);
    int32_t hgid = bam_name2id(hdr, hreg.chr.c_str()); // hgene chr idx
    int32_t tgid = bam_name2id(hdr, treg.chr.c_str()); // tgene chr idx
    int32_t hpos = std::atoi(hreg.beg.c_str());
    int32_t hend = std::atoi(hreg.end.c_str());
    int32_t tpos = std::atoi(treg.beg.c_str());
    int32_t tend = std::atoi(treg.end.c_str());
    int32_t srcnt = 0, oscnt = 0, dpcnt = 0;
    samFile* srfp = sam_open(srbam.c_str(), "wb");
    samFile* osfp = NULL;
    if(!osbam.empty()){
        sam_open(osbam.c_str(), "wb");
        assert(sam_hdr_write(osfp, hdr) >= 0);
    }
    samFile* dpfp = sam_open(dpbam.c_str(), "wb");
    assert(sam_hdr_write(srfp, hdr) >= 0 && sam_hdr_write(dpfp, hdr) >= 0);
    hts_idx_t* idx = sam_index_load(sfp, inbam.c_str());
    // check hgene range for possible evidence
     const uint16_t BAM_SRSKIP_MASK = (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY);
    hts_itr_t* itr = sam_itr_querys(idx, hdr, hqreg.c_str());
    bam1_t* b = bam_init1();
    while(sam_itr_next(sfp, itr, b) >= 0){
        if(b->core.flag & BAM_SRSKIP_MASK) continue;
        std::pair<int32_t, int32_t> clips = bamutil::getSoftClipLength(b);
        if((clips.first  > 0) ^ (clips.second > 0)){
            uint8_t* sa = bam_aux_get(b, "SA");
            if(sa){
                std::string sastr = bam_aux2Z(sa);
                std::vector<std::string> cvs;
                std::vector<std::string> vstr;
                util::split(sastr, cvs, ";");
                if(cvs[1].empty()) sastr = cvs[0];
                else{
                    std::vector<int32_t> mvidx;
                    for(uint32_t cvidx = 0; cvidx < cvs.size() - 1; ++cvidx){
                        util::split(cvs[cvidx], vstr, ",");
                        if(vstr[3].find_first_of("SH") == vstr[3].find_last_of("SH")){
                            mvidx.push_back(cvidx);
                        }
                    }
                    if(mvidx.size() == 1){
                        sastr = cvs[mvidx[0]];
                    }else{
                        sastr = "";
                    }
                }
                if(sastr.empty()) continue;
                util::split(sastr, vstr, ",");
                std::string sachr = vstr[0];
                int32_t sapos = std::atoi(vstr[1].c_str());
                if(sachr == treg.chr && sapos > tpos && sapos < tend){
                    ++srcnt;
                    assert(sam_write1(srfp, hdr, b) >= 0);
                }
            }else{
                if(osfp){
                    ++oscnt;
                    assert(sam_write1(osfp, hdr, b) >= 0);
                }
            }
        }
        if(b->core.mtid == tgid){
            if(b->core.mpos > tpos && b->core.mpos < tend){
                ++dpcnt;
                assert(sam_write1(dpfp, hdr, b) >= 0);
            }
        }
    }
    // check tgene range for possible evidence
    hts_itr_destroy(itr);
    itr = sam_itr_querys(idx, hdr, tqreg.c_str());
    while(sam_itr_next(sfp, itr, b) >= 0){
        if(b->core.flag & BAM_SRSKIP_MASK) continue;
        std::pair<int32_t, int32_t> clips = bamutil::getSoftClipLength(b);
        if((clips.first > 0) ^ (clips.second > 0)){
            uint8_t* sa = bam_aux_get(b, "SA");
            if(sa){
                std::string sastr = bam_aux2Z(sa);
                std::vector<std::string> cvs;
                std::vector<std::string> vstr;
                util::split(sastr, cvs, ";");
                if(cvs[1].empty()) sastr = cvs[0];
                else{
                    std::vector<int32_t> mvidx;
                    for(uint32_t cvidx = 0; cvidx < cvs.size() - 1; ++cvidx){
                        util::split(cvs[cvidx], vstr, ",");
                        if(vstr[3].find_first_of("SH") == vstr[3].find_last_of("SH")){
                            mvidx.push_back(cvidx);
                        }
                    }
                    if(mvidx.size() == 1){
                        sastr = cvs[mvidx[0]];
                    }else{
                        sastr = "";
                    }
                }
                if(sastr.empty()) continue;
                util::split(sastr, vstr, ",");
                std::string sachr = vstr[0];
                int32_t sapos = std::atoi(vstr[1].c_str());
                if(sachr == hreg.chr && sapos > hpos && sapos < hend){
                    ++srcnt;
                    assert(sam_write1(srfp, hdr, b) >= 0);
                }
            }else{
                if(osfp){
                    ++oscnt;
                    assert(sam_write1(osfp, hdr, b) >= 0);
                }
            }
        }
        if(b->core.mtid == hgid){
            if(b->core.mpos > hpos && b->core.mpos < hend){
                ++dpcnt;
                assert(sam_write1(dpfp, hdr, b) >= 0);
            }
        }
    }
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    sam_close(sfp);
    sam_close(srfp);
    if(osfp) sam_close(osfp);
    sam_close(dpfp);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
    // write statics to output 
    std::ofstream fw(table);
    fw << "fusion" << "\t" << hgene << "->" << tgene << "\n";
    fw << "srcount" << "\t" << srcnt << "\n";
    if(oscnt) fw << "oscount" << "\t" << oscnt << "\n";
    fw << "dpcount" << "\t" << dpcnt << "\n";
    fw.close();
}

void SVDebug::debugRNA(){
    // fetch gene coordinates
    std::map<std::string, std::string> g2tmap;
    util::makeMapPairFromFileByLine(gene2trans, g2tmap);
    std::string htrs = g2tmap[hgene];
    std::string ttrs = g2tmap[tgene];
    // open bam and idx
    samFile* sfp = sam_open(inbam.c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(sfp);
    int32_t hgid = bam_name2id(hdr, htrs.c_str()); // hgene chr idx
    int32_t tgid = bam_name2id(hdr, ttrs.c_str()); // tgene chr idx
    int32_t hpos = 0;
    int32_t tpos = 0;
    int32_t hend = hdr->target_len[hgid];
    int32_t tend = hdr->target_len[tgid];
    int32_t srcnt = 0, oscnt = 0, dpcnt = 0;
    samFile* srfp = sam_open(srbam.c_str(), "wb");
    samFile* osfp = NULL;
    if(!osbam.empty()){
        sam_open(osbam.c_str(), "wb");
        assert(sam_hdr_write(osfp, hdr) >= 0);
    }
    samFile* dpfp = sam_open(dpbam.c_str(), "wb");
    assert(sam_hdr_write(srfp, hdr) >= 0 && sam_hdr_write(dpfp, hdr) >= 0);
    hts_idx_t* idx = sam_index_load(sfp, inbam.c_str());
    // check hgene range for possible evidence
    const uint16_t BAM_SRSKIP_MASK = (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY);
    std::stringstream regs;
    regs << htrs << ":" << hpos << "-" << hend;
    hts_itr_t* itr = sam_itr_querys(idx, hdr, regs.str().c_str());
    bam1_t* b = bam_init1();
    while(sam_itr_next(sfp, itr, b) >= 0){
        if(b->core.flag & BAM_SRSKIP_MASK) continue;
        std::pair<int32_t, int32_t> clips = bamutil::getSoftClipLength(b);
        if(clips.first || clips.second){
            uint8_t* data = bam_aux_get(b, "SA");
            if(data){
                char* val = bam_aux2Z(data);
                std::vector<std::string> vstr;
                util::split(val, vstr, ",");
                std::string sachr = vstr[0];
                int32_t sapos = std::atoi(vstr[1].c_str());
                if(sachr == ttrs && sapos > tpos && sapos < tend){
                    ++srcnt;
                    assert(sam_write1(srfp, hdr, b) >= 0);
                }
            }else{
                if(osfp){
                    ++oscnt;
                    assert(sam_write1(osfp, hdr, b) >= 0);
                }
            }
        }
        if(b->core.mtid == tgid){
            if(b->core.mpos > tpos && b->core.mpos < tend){
                ++dpcnt;
                assert(sam_write1(dpfp, hdr, b) >= 0);
            }
        }
    }
    // check tgene range for possible evidence
    hts_itr_destroy(itr);
    regs.clear();
    regs.str("");
    regs << ttrs << ":" << tpos << "-" << tend;
    itr = sam_itr_querys(idx, hdr, regs.str().c_str());
    while(sam_itr_next(sfp, itr, b) >= 0){
        if(b->core.flag & BAM_SRSKIP_MASK) continue;
        std::pair<int32_t, int32_t> clips = bamutil::getSoftClipLength(b);
        if(clips.first || clips.second){
            uint8_t* data = bam_aux_get(b, "SA");
            if(data){
                char* val = bam_aux2Z(data);
                std::vector<std::string> vstr;
                util::split(val, vstr, ",");
                std::string sachr = vstr[0];
                int32_t sapos = std::atoi(vstr[1].c_str());
                if(sachr == htrs && sapos > hpos && sapos < hend){
                    ++srcnt;
                    assert(sam_write1(srfp, hdr, b) >= 0);
                }
            }else{
                if(osfp){
                    ++oscnt;
                    assert(sam_write1(osfp, hdr, b) >= 0);
                }
            }
        }
        if(b->core.mtid == hgid){
            if(b->core.mpos > hpos && b->core.mpos < hend){
                ++dpcnt;
                assert(sam_write1(dpfp, hdr, b) >= 0);
            }
        }
    }
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    sam_close(sfp);
    sam_close(srfp);
    if(osfp) sam_close(osfp);
    sam_close(dpfp);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
    // write statics to output 
    std::ofstream fw(table);
    fw << "fusion" << "\t" << hgene << "->" << tgene << "\n";
    fw << "srcount" << "\t" << srcnt << "\n";
    if(oscnt) fw << "oscount" << "\t" << oscnt << "\n";
    fw << "dpcount" << "\t" << dpcnt << "\n";
    fw.close();
}

