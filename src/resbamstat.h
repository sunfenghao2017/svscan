#ifndef RES_BAM_STAT_H
#define RES_BAM_STAT_H

#include "realnfilter.h"
#include <htslib/sam.h>
#include <cstdint>
#include <string>
#include <map>

/** structure to store pe partner support status */
struct PePtnStat{
    bool is_read1 = false;
    int32_t svid = -1;
    bool found = false;
    bool skip = false;
    int mq = 0;
    std::string chr;
    int32_t mpos = -1;
    bool valid = false;
};

typedef std::map<std::string, PePtnStat*> PePtnMap;

// fusion reads pattern count index map
const std::map<std::string, int> FsPatMatIdx = {{"++", 0}, {"+-", 1}, {"--", 2}, {"-+", 3}};

/** class to store a read supporting status */
struct ReadSupport{
    int32_t mR1SVID = -1;
    uint8_t mR1MapQ = 0;
    int8_t mR1SRT = -1;
    int32_t mR1Hit = 0;
    int32_t mR1PHit = 0;
    int32_t mR1PSPos = -1;
    int32_t mR1PTid = -1;
    std::string mR1PChr;
    std::string mR1PStrand = "-";
    int32_t mR1PBeg = -1;
    int32_t mR1PEnd = -1;
    int mR1PTgt = 0;
    int32_t mR1SHit = 0;
    int32_t mR1SSPos = -1;
    int32_t mR1STid = -1;
    std::string mR1SChr;
    std::string mR1SStrand = "-";
    int32_t mR1SBeg = -1;
    int32_t mR1SEnd = -1;
    int mR1STgt = 0;
    int mR1Seed = 0;
    
    int32_t mR2SVID = -1;
    uint8_t mR2MapQ = 0;
    int8_t mR2SRT = -1;
    int32_t mR2Hit = 0;
    int32_t mR2PHit = 0;
    int32_t mR2PSPos = -1;
    int32_t mR2PTid = -1;
    std::string mR2PChr;
    std::string mR2PStrand = "-";
    int32_t mR2PBeg = -1;
    int32_t mR2PEnd = -1;
    int mR2PTgt = 0;
    int32_t mR2SHit = 0;
    int32_t mR2SSPos = -1;
    int32_t mR2STid = -1;
    std::string mR2SChr;
    std::string mR2SStrand = "-";
    int32_t mR2SBeg = -1;
    int32_t mR2SEnd = -1;
    int mR2STgt = 0;
    int mR2Seed = 0;

    inline friend std::ostream& operator<<(std::ostream& os, const ReadSupport& rs){
        os << "==========read1==========\n";
        os << "mR1SVID: " << rs.mR1SVID << "\n";
        os << "mR1MapQ: " << rs.mR1MapQ << "\n";
        os << "mR1SRT:  " << rs.mR1SRT  << "\n";
        os << "mR1Hit:  " << rs.mR1Hit  << "\n";
        os << "mR1PHit: " << rs.mR1PHit << "\n";
        os << "mR1PTid: " << rs.mR1PTid << "\n";
        os << "mR1PChr: " << rs.mR1PChr << "\n";
        os << "mR1PSPos:" << rs.mR1PSPos << "\n";
        os << "mR1PStrd:" << rs.mR1PStrand << "\n";
        os << "mR1PBeg: " << rs.mR1PBeg << "\n";
        os << "mR1PEnd: " << rs.mR1PEnd << "\n";
        os << "mR1PTgt: " << rs.mR1PTgt << "\n";
        os << "mR1SHit: " << rs.mR1SHit << "\n";
        os << "mR1STid: " << rs.mR1STid << "\n";
        os << "mR1SChr: " << rs.mR1SChr << "\n";
        os << "mR1SSPos:" << rs.mR1SSPos << "\n";
        os << "mR1SStrd:" << rs.mR1SStrand << "\n";
        os << "mR1SBeg: " << rs.mR1SBeg << "\n";
        os << "mR1SEnd: " << rs.mR1SEnd << "\n";
        os << "mR1STgt: " << rs.mR1STgt << "\n";
        os << "mR1Seed: " << rs.mR1Seed << "\n";
        os << "==========read2==========\n";
        os << "mR2SVID: " << rs.mR2SVID << "\n";
        os << "mR2MapQ: " << rs.mR2MapQ << "\n";
        os << "mR2SRT:  " << rs.mR2SRT  << "\n";
        os << "mR2Hit:  " << rs.mR2Hit  << "\n";
        os << "mR2PHit: " << rs.mR2PHit << "\n";
        os << "mR2PTid: " << rs.mR2PTid << "\n";
        os << "mR2PChr: " << rs.mR2PChr << "\n";
        os << "mR2PSPos:" << rs.mR2PSPos << "\n";
        os << "mR2PStrd:" << rs.mR2PStrand << "\n";
        os << "mR2PBeg: " << rs.mR2PBeg << "\n";
        os << "mR2PEnd: " << rs.mR2PEnd << "\n";
        os << "mR2PTgt: " << rs.mR2PTgt << "\n";
        os << "mR2SHit: " << rs.mR2SHit << "\n";
        os << "mR2STid: " << rs.mR2STid << "\n";
        os << "mR2SChr: " << rs.mR2SChr << "\n";
        os << "mR2SSPos:" << rs.mR2SSPos << "\n";
        os << "mR2SStrd:" << rs.mR2SStrand << "\n";
        os << "mR2SBeg: " << rs.mR2SBeg << "\n";
        os << "mR2SEnd: " << rs.mR2SEnd << "\n";
        os << "mR2STgt: " << rs.mR2STgt << "\n";
        os << "mR2Seed: " << rs.mR2Seed << "\n";
        return os;
    }

    inline void countPattern(int32_t chr1, int32_t pos1, int& r1p, int& r2p, bool issrsv){
        r1p = -1; r2p = -1; // initialize to invalid pattern
        bool phitpos1 = false;
        if(mR1Seed == 1){
            if(mR1PTid == chr1){
                if(mR1STid != chr1){
                    phitpos1 = true;
                }else{
                    if(std::abs(mR1PSPos - pos1) < std::abs(mR1SSPos - pos1)) phitpos1 = true;
                }
            }
            std::string pat;
            if(phitpos1) pat = mR1PStrand + mR1SStrand;
            else pat = mR1SStrand + mR1PStrand; 
            auto iter = FsPatMatIdx.find(pat);
            r1p = iter->second;
        }else if(mR1SRT == 1 && (!issrsv)){
            if(mR1PTid == chr1){
                if(mR2PTid != chr1){
                    phitpos1 = true;
                }else{
                    if(std::abs(mR1PSPos - pos1) < std::abs(mR2PSPos - pos1)) phitpos1 = true;
                }
            }
            std::string pat;
            if(phitpos1) pat = mR1PStrand + mR2PStrand;
            else pat = mR2PStrand + mR1PStrand;
            auto iter = FsPatMatIdx.find(pat);
            r1p = iter->second;
        }

        phitpos1 = false;
        if(mR2Seed == 1){
            if(mR2PTid == chr1){
                if(mR2STid != chr1){
                    phitpos1 = true;
                }else{
                    if(std::abs(mR2PSPos - pos1) < std::abs(mR2SSPos - pos1)) phitpos1 = true;
                }
            }
            std::string pat;
            if(phitpos1) pat = mR2PStrand + mR2SStrand;
            else pat = mR2SStrand + mR2PStrand;
            auto iter = FsPatMatIdx.find(pat);
            r2p = iter->second;
        }else if(mR2SRT == 1 && (mR1SRT != 1) && (!issrsv)){
            if(mR2PTid != chr1){
                phitpos1 = true;
            }else{
                if(std::abs(mR1PSPos - pos1) < std::abs(mR2PSPos - pos1)) phitpos1 = true;
            }
            std::string pat;
            if(phitpos1) pat = mR1PStrand + mR2PStrand;
            else pat = mR2PStrand + mR1PStrand;
            auto iter = FsPatMatIdx.find(pat);
            r2p = iter->second;
        }
    }

};

/** type to store read supporting statistics */
typedef std::map<std::string, ReadSupport*> ReadSupportStatMap;

inline void getReadSupportStatus(const std::string& bam, ReadSupportStatMap& rssm, PePtnMap& pem,  RealnFilter* rf){
    samFile* fp = sam_open(bam.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    while(sam_read1(fp, h, b) >= 0){
        uint8_t* idata = bam_aux_get(b, "ZF");
        uint8_t* sdata = bam_aux_get(b, "ST");
        uint8_t* saval = bam_aux_get(b, "SA");
        int32_t phit = 0, pbeg = -1, pend = -1, shit = 0, sbeg = -1, send = -1, pspos = -1, sspos = -1;
        int32_t ptid = -1, stid = -1;
        std::string pseq, sseq;
        std::string pchr, schr;
        std::string pstd, sstd;
        if(saval){
            ptid = b->core.tid;
            if(b->core.flag & BAM_FREVERSE) pstd = "-";
            else pstd = "+";
            // get optimal SA
            std::string sastr = bam_aux2Z(saval);
            std::vector<std::string> cvs;
            std::vector<std::string> vstr;
            util::split(sastr, cvs, ";");
            std::string optsa;
            if(cvs[1].empty()){
                util::split(cvs[0], vstr, ",");
                schr = vstr[0];
                stid = sam_hdr_name2tid(h, schr.c_str());
                sbeg = std::atoi(vstr[1].c_str()) - 1;
                sspos = sbeg;
                send = sbeg;
                sstd = vstr[2];
                optsa = vstr[3];
            }else{
                for(uint32_t cvidx = 0; cvidx < cvs.size() - 1; ++cvidx){
                    util::split(cvs[cvidx], vstr, ",");
                    if(vstr[3].find_first_of("SH") == vstr[3].find_last_of("SH")){
                        schr = vstr[0];
                        stid = sam_hdr_name2tid(h, schr.c_str());
                        sbeg = std::atoi(vstr[1].c_str()) - 1;
                        sspos = sbeg;
                        send = sbeg;
                        sstd = vstr[2];
                        optsa = vstr[3];
                        break;
                    }
                }
            }
            // parse primary alignment length
            int32_t sclen = 0;
            bool leadsc = false;
            for(uint32_t i = 0; i < b->core.n_cigar; ++i){
                if(bam_cigar_op(bam_get_cigar(b)[i]) == BAM_CSOFT_CLIP){
                    sclen = bam_cigar_oplen(bam_get_cigar(b)[i]);
                    if(i == 0) leadsc = true;
                    break;
                }
            }
            if(leadsc){
                pspos = b->core.pos;
                pend = bam_endpos(b) - 1;
            }else{
                pspos = bam_endpos(b) - 1;
                pend = pspos;
            }
            // parse supplementary alignment length
            int32_t sascl = 0;
            char* scg = const_cast<char*>(optsa.c_str());
            while(*scg && *scg != '*'){
                if(isdigit((int)*scg)){
                    long num = strtol(scg, &scg, 10);
                    if(*scg == 'S'){
                        sascl = num;
                        if(send != sbeg){
                            --send;
                            --sspos;
                        }
                        break;
                    }else if(*scg == 'M' || *scg == 'X' || *scg == '=' || *scg == 'D'){
                        sspos += num;
                        send += num;
                    }
                }
                ++scg;
            }
            // get seq
            pchr = h->target_name[b->core.tid];
            pbeg = b->core.pos;
            std::string rseq(b->core.l_qseq, '\0');
            for(int32_t i = 0; i < b->core.l_qseq; ++i){
                rseq[i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];
            }
            if(leadsc){
                pseq = rseq.substr(sclen);
                sseq = rseq.substr(0, sclen);
            }else{
                pseq = rseq.substr(0, b->core.l_qseq - sclen);
                sseq = rseq.substr(b->core.l_qseq - sclen);
            }
            bool pfm = false, sfm = false;
            int32_t pftid = -1, pfb = -1, pfe = -1, sftid = -1, sfb = -1, sfe = -1;
            phit = rf->validSRSeq(pseq, pfm, pftid, pfb, pfe);
            shit = rf->validSRSeq(sseq, sfm, sftid, sfb, sfe);
            int32_t maxoff = std::max(10, b->core.l_qseq - sclen - sascl);
            if(phit == 1 && pfm){// check fm
                if(b->core.tid == pftid &&
                   std::abs(b->core.pos - pfb) < maxoff &&
                   std::abs(pend - pfe) < maxoff){
                    // donothing
                }else{
                    phit = 1000;
                }
            }
            if(shit == 1 && sfm){// check fm
                if(stid == sftid &&
                   std::abs(sbeg - sfb) < maxoff &&
                   std::abs(send - sfe) < maxoff){
                    // donothing
                }else{
                    shit = 1000;
                }
            }
            if(shit == 1){// check ins removed part
                int inslen = sseq.length() - (b->core.l_qseq - sascl);
                if(inslen > 0){
                    shit = rf->validSRSeq(sseq.substr(inslen), sfm, sftid, sfb, sfe);
                }
            }
        }
        if(idata && sdata){
            int svid = bam_aux2i(idata);
            int srst = bam_aux2i(sdata);
            std::string qname = bam_get_qname(b);
            auto iter = rssm.find(qname);
            if(iter == rssm.end()){
                ReadSupport* rs = new ReadSupport();
                if(b->core.flag & BAM_FREAD1){
                    rs->mR1SVID = svid;
                    rs->mR1MapQ = b->core.qual;
                    rs->mR1SRT = srst;
                    rs->mR1Hit = phit + shit;
                    if(saval){
                        rs->mR1PHit = phit;
                        rs->mR1SHit = shit;
                        rs->mR1PSPos = pspos;
                        rs->mR1SSPos = sspos;
                        rs->mR1PTid = ptid;
                        rs->mR1STid = stid;
                        rs->mR1PChr = pchr;
                        rs->mR1SChr = schr;
                        rs->mR1PStrand = pstd;
                        rs->mR1SStrand = sstd;
                        rs->mR1PBeg = pbeg;
                        rs->mR1PEnd = pend;
                        rs->mR1SBeg = sbeg;
                        rs->mR1SEnd = send;
                        rs->mR1Seed = 1;
                    }else if(srst){ // dp type
                        rs->mR1PTid = b->core.tid;
                        rs->mR2PTid = b->core.mtid;
                        rs->mR1PSPos = b->core.pos;
                        rs->mR2PSPos = b->core.mpos;
                        if(b->core.flag & BAM_FREVERSE){
                            rs->mR1PStrand = "-";
                        }else{
                            rs->mR1PStrand = "+";
                        }
                        if(b->core.flag & BAM_FMREVERSE){
                            rs->mR2PStrand = "-";
                        }else{
                            rs->mR2PStrand = "+";
                        }
                    }
                }else{
                    rs->mR2SVID = svid;
                    rs->mR2MapQ = b->core.qual;
                    rs->mR2SRT = srst;
                    rs->mR2Hit = phit + shit;
                    if(saval){
                        rs->mR2PHit = phit;
                        rs->mR2SHit = shit;
                        rs->mR2PSPos = pspos;
                        rs->mR2SSPos = sspos;
                        rs->mR2PTid = ptid;
                        rs->mR2STid = stid;
                        rs->mR2PChr = pchr;
                        rs->mR2SChr = schr;
                        rs->mR2PStrand = pstd;
                        rs->mR2SStrand = sstd;
                        rs->mR2PBeg = pbeg;
                        rs->mR2PEnd = pend;
                        rs->mR2SBeg = sbeg;
                        rs->mR2SEnd = send;
                        rs->mR2Seed = 1;
                    }else if(srst){
                        rs->mR2PTid = b->core.tid;
                        rs->mR1PTid = b->core.mtid;
                        rs->mR2PSPos = b->core.pos;
                        rs->mR1PSPos = b->core.mpos;
                        if(b->core.flag & BAM_FREVERSE){
                            rs->mR2PStrand = "-";
                        }else{
                            rs->mR2PStrand = "+";
                        }
                        if(b->core.flag & BAM_FMREVERSE){
                            rs->mR1PStrand = "-";
                        }else{
                            rs->mR1PStrand = "+";
                        }
                    }
                }
                rssm[qname] = rs;
            }else{
                if(b->core.flag & BAM_FREAD1){
                    iter->second->mR1SVID = svid;
                    iter->second->mR1MapQ = b->core.qual;
                    iter->second->mR1SRT = srst;
                    iter->second->mR1Hit = phit + shit;
                    if(saval){
                        iter->second->mR1PHit = phit;
                        iter->second->mR1SHit = shit;
                        iter->second->mR1PSPos = pspos;
                        iter->second->mR1SSPos = sspos;
                        iter->second->mR1PTid = ptid;
                        iter->second->mR1STid = stid;
                        iter->second->mR1PChr = pchr;
                        iter->second->mR1SChr = schr;
                        iter->second->mR1PStrand = pstd;
                        iter->second->mR1SStrand = sstd;
                        iter->second->mR1PBeg = pbeg;
                        iter->second->mR1PEnd = pend;
                        iter->second->mR1SBeg = sbeg;
                        iter->second->mR1SEnd = send;
                        iter->second->mR1Seed = 1;
                    }else if(srst){
                        iter->second->mR1PTid = b->core.tid;
                        iter->second->mR2PTid = b->core.mtid;
                        iter->second->mR1PSPos = b->core.pos;
                        iter->second->mR2PSPos = b->core.mpos;
                        if(b->core.flag & BAM_FREVERSE){
                            iter->second->mR1PStrand = "-";
                        }else{
                            iter->second->mR1PStrand = "+";
                        }
                        if(b->core.flag & BAM_FMREVERSE){
                            iter->second->mR2PStrand = "-";
                        }else{
                            iter->second->mR2PStrand = "+";
                        }
                    }
                }else{
                    iter->second->mR2SVID = svid;
                    iter->second->mR2MapQ = b->core.qual;
                    iter->second->mR2SRT = srst;
                    iter->second->mR2Hit = phit + shit;
                    if(saval){
                        iter->second->mR2PHit = phit;
                        iter->second->mR2SHit = shit;
                        iter->second->mR2PSPos = pspos;
                        iter->second->mR2SSPos = sspos;
                        iter->second->mR2PTid = ptid;
                        iter->second->mR2STid = stid;
                        iter->second->mR2PChr = pchr;
                        iter->second->mR2SChr = schr;
                        iter->second->mR2PStrand = pstd;
                        iter->second->mR2SStrand = sstd;
                        iter->second->mR2PBeg = pbeg;
                        iter->second->mR2PEnd = pend;
                        iter->second->mR2SBeg = sbeg;
                        iter->second->mR2SEnd = send;
                        iter->second->mR2Seed = 1;
                    }else if(srst){
                        iter->second->mR2PTid = b->core.tid;
                        iter->second->mR1PTid = b->core.mtid;
                        iter->second->mR2PSPos = b->core.pos;
                        iter->second->mR1PSPos = b->core.mpos;
                        if(b->core.flag & BAM_FREVERSE){
                            iter->second->mR2PStrand = "-";
                        }else{
                            iter->second->mR2PStrand = "+";
                        }
                        if(b->core.flag & BAM_FMREVERSE){
                            iter->second->mR1PStrand = "-";
                        }else{
                            iter->second->mR1PStrand = "+";
                        }
                    }
                }
            }
            if(srst == 1){//collect pe
                PePtnStat *pps = new PePtnStat();
                pps->found = false;
                pps->skip = false;
                pps->svid = svid;
                if(b->core.flag & BAM_FREAD1) pps->is_read1 = false;
                else pps->is_read1 = true;
                pps->chr = sam_hdr_tid2name(h, b->core.mtid);
                pps->mpos = b->core.mpos;
                pps->mq = b->core.qual;
                pps->valid = false;
                pem[bam_get_qname(b)] = pps;
            }
        }
    }
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
}

#endif
