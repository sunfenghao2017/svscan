#include "rsvg.h"

void RSVGen::getunit(SeqInfo& si, std::vector<UnitRec>& vur){
    hts_itr_t* itr = tbx_itr_querys(tbx, si.trs.c_str());
    kstring_t rec = {0, 0, 0};
    std::vector<std::string> vstr;
    while(tbx_itr_next(tfp, tbx, itr, &rec) >= 0){
        util::split(rec.s, vstr, "\t");
        if(util::startsWith(vstr[3], "utr")) continue;
        if(si.chr.empty()) si.chr = vstr[6];
        UnitRec ur;
        ur.rbeg = std::atoi(vstr[1].c_str());
        ur.rend = std::atoi(vstr[2].c_str());
        ur.gbeg = std::atoi(vstr[7].c_str());
        ur.gend = std::atoi(vstr[8].c_str());
        ur.len = ur.rend - ur.rbeg + 1;
        ur.cnt = std::atoi(vstr[4].c_str());
        vur.push_back(ur);
    }
    // cleanup
    hts_itr_destroy(itr);
    ks_release(&rec);
}

void RSVGen::gethseq(SeqInfo& si){
    // fetch all exon of this transcript
    std::vector<UnitRec> vur;
    getunit(si, vur);
    // compute left flank exon surround this transcript
    si.dend = vur[si.exon].gend;
    si.rend = vur[si.exon].rend;
    si.dbeg = vur[si.exon].gbeg;
    si.rbeg = vur[si.exon].rbeg;
    int i = si.exon - 1;
    while(i >= 0){
        if(si.flklen < maxflk){
            if(si.flklen + vur[i].len <= maxflk){
                si.flklen += vur[i].len;
                si.dbeg = vur[i].gbeg;
                si.rbeg = vur[i].rbeg;
            }else{
                int offbp = si.flklen + vur[i].len - maxflk;
                si.flklen = maxflk;
                si.dbeg = vur[i].gbeg + offbp;
                si.rbeg = vur[i].rbeg + offbp;
            }
        }
        --i;
    }
    // fetch sequence
    int sl = -1;
    si.seq = faidx_fetch_seq(fai, si.trs.c_str(), si.rbeg, si.rend, &sl);
}

void RSVGen::gettseq(SeqInfo& si){
    // fetch all exon of this transcript
    std::vector<UnitRec> vur;
    getunit(si, vur);
    // compute left flank exon surround this transcript
    si.dend = vur[si.exon].gend;
    si.rend = vur[si.exon].rend;
    si.dbeg = vur[si.exon].gbeg;
    si.rbeg = vur[si.exon].rbeg;
    int i = si.exon + 1;
    int m = vur.size();
    while(i < m){
        if(si.flklen < maxflk){
            if(si.flklen + vur[i].len <= maxflk){
                si.flklen += vur[i].len;
                si.dend = vur[i].gend;
                si.rend = vur[i].rend;
            }else{
                int offbp = si.flklen + vur[i].len - maxflk;
                si.flklen = maxflk;
                si.dend = vur[i].gend - offbp;
                si.rend = vur[i].rend - offbp;
            }
        }
        ++i;
    }
    // fetch sequence
    int sl = -1;
    si.seq = faidx_fetch_seq(fai, si.trs.c_str(), si.rbeg, si.rend, &sl);
}

void RSVGen::gensv(){
    // parse input config file
    std::ifstream fr(incfg);
    std::vector<std::string> vrec;
    std::vector<FuseInfo> vfi;
    std::string line;
    while(std::getline(fr, line)){
        util::split(line, vrec, "\t");
        // FuseInfo
        FuseInfo fi;
        // hseq
        fi.hseq.gene = vrec[0];
        fi.hseq.trs = g2tm[fi.hseq.gene];
        fi.hseq.exon = std::atoi(vrec[1].c_str()) - 1;
        // tseq
        fi.tseq.gene = vrec[2];
        fi.tseq.trs = g2tm[fi.tseq.gene];
        fi.tseq.exon = std::atoi(vrec[3].c_str()) - 1;
        // push
        vfi.push_back(fi);
    }
    fr.close();
    // get sequence
    for(uint32_t i = 0; i < vfi.size(); ++i){
        gethseq(vfi[i].hseq);
        gettseq(vfi[i].tseq);
    }
    // output result
    std::ofstream fwfa(outfa), fwcf(outcfg);
    FuseInfo::outInfoHead(fwcf);
    for(uint32_t i = 0; i < vfi.size(); ++i){
        vfi[i].outSeq(fwfa);
        vfi[i].outInfoRec(fwcf);
    }
    fwfa.close();
    fwcf.close();
}
