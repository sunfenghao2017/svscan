#ifndef RES_BAM_STAT_H
#define RES_BAM_STAT_H

#include "realnfilter.h"
#include <htslib/sam.h>
#include <cstdint>
#include <string>
#include <map>

/** class to store a read supporting status */
struct ReadSupport{
    int32_t mR1SVID = -1;
    uint8_t mR1MapQ = 0;
    int8_t mR1SRT = -1;
    int32_t mR1Hit = 0;
    int32_t mR1PHit = 0;
    std::string mR1PChr;
    int32_t mR1PBeg = -1;
    int32_t mR1PEnd = -1;
    int mR1PTgt = 0;
    int32_t mR1SHit = 0;
    std::string mR1SChr;
    int32_t mR1SBeg = -1;
    int32_t mR1SEnd = -1;
    int mR1STgt = 0;
    int mR1Seed = 0;
    
    int32_t mR2SVID = -1;
    uint8_t mR2MapQ = 0;
    int8_t mR2SRT = -1;
    int32_t mR2Hit = 0;
    int32_t mR2PHit = 0;
    std::string mR2PChr;
    int32_t mR2PBeg = -1;
    int32_t mR2PEnd = -1;
    int mR2PTgt = 0;
    int32_t mR2SHit = 0;
    std::string mR2SChr;
    int32_t mR2SBeg = -1;
    int32_t mR2SEnd = -1;
    int mR2STgt = 0;
    int mR2Seed = 0;

    inline friend std::ostream& operator<<(std::ostream& os, const ReadSupport& rs){
        os << "==========read1==========\n";
        os << "mR1SVID: " << rs.mR1SVID << "\n";
        os << "mR1MapQ: " << rs.mR1MapQ << "\n";
        os << "mR1SRT:  " << rs.mR1SRT << "\n";
        os << "mR1Hit:  " << rs.mR1Hit << "\n";
        os << "mR1PHit: " << rs.mR1PHit << "\n";
        os << "mR1PChr: " << rs.mR1PChr << "\n";
        os << "mR1PBeg: " << rs.mR1PBeg << "\n";
        os << "mR1PEnd: " << rs.mR1PEnd << "\n";
        os << "mR1PTgt: " << rs.mR1PTgt << "\n";
        os << "mR1SHit: " << rs.mR1SHit << "\n";
        os << "mR1SChr: " << rs.mR1SChr << "\n";
        os << "mR1SBeg: " << rs.mR1SBeg << "\n";
        os << "mR1SEnd: " << rs.mR1SEnd << "\n";
        os << "mR1STgt: " << rs.mR1STgt << "\n";
        os << "mR1Seed: " << rs.mR1Seed << "\n";
        os << "==========read2==========\n";
        os << "mR2SVID: " << rs.mR2SVID << "\n";
        os << "mR2MapQ: " << rs.mR2MapQ << "\n";
        os << "mR2SRT:  " << rs.mR2SRT << "\n";
        os << "mR2Hit:  " << rs.mR2Hit << "\n";
        os << "mR2PHit: " << rs.mR2PHit << "\n";
        os << "mR2PChr: " << rs.mR2PChr << "\n";
        os << "mR2PBeg: " << rs.mR2PBeg << "\n";
        os << "mR2PEnd: " << rs.mR2PEnd << "\n";
        os << "mR2PTgt: " << rs.mR2PTgt << "\n";
        os << "mR2SHit: " << rs.mR2SHit << "\n";
        os << "mR2SChr: " << rs.mR2SChr << "\n";
        os << "mR2SBeg: " << rs.mR2SBeg << "\n";
        os << "mR2SEnd: " << rs.mR2SEnd << "\n";
        os << "mR2STgt: " << rs.mR2STgt << "\n";
        os << "mR2Seed: " << rs.mR2Seed << "\n";
        return os;
    }
};

/** type to store read supporting statistics */
typedef std::map<std::string, ReadSupport*> ReadSupportStatMap;

inline void getReadSupportStatus(const std::string& bam, ReadSupportStatMap& rssm, RealnFilter* rf){
    samFile* fp = sam_open(bam.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    while(sam_read1(fp, h, b) >= 0){
        uint8_t* idata = bam_aux_get(b, "ZF");
        uint8_t* sdata = bam_aux_get(b, "ST");
        uint8_t* saval = bam_aux_get(b, "SA");
        int32_t phit = 0, pbeg = -1, pend = -1, shit = 0, sbeg = -1, send = -1;
        std::string pchr, schr;
        if(saval){
            int32_t sclen = 0;
            bool leadsc = false;
            for(uint32_t i = 0; i < b->core.n_cigar; ++i){
                if(bam_cigar_op(bam_get_cigar(b)[i]) == BAM_CSOFT_CLIP){
                    sclen = bam_cigar_oplen(bam_get_cigar(b)[i]);
                    if(i == 0) leadsc = true;
                    break;
                }
            }
            pchr = h->target_name[b->core.tid];
            pbeg = b->core.pos;
            pend = bam_endpos(b);
            std::string qseq(b->core.l_qseq, '\0');
            for(int32_t i = 0; i < b->core.l_qseq; ++i){
                qseq[i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];
            }
            std::string pseq, sseq;
            if(leadsc){
                sseq = qseq.substr(0, sclen);
                pseq = qseq.substr(sclen);
            }else{
                sseq = qseq.substr(sclen);
                pseq = qseq.substr(0, sclen);
            }
            phit = rf->validSRSeq(pseq);
            shit = rf->validSRSeq(sseq);
            std::string ss = bam_aux2Z(saval);
            std::vector<std::string> vsas;
            util::split(ss, vsas, ";");
            std::vector<std::string> vstr;
            util::split(vsas[0], vstr, ",");
            schr = vstr[0];
            sbeg = std::atoi(vstr[1].c_str());
            send = sbeg + sclen;
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
                        rs->mR1PChr = pchr;
                        rs->mR1SChr = schr;
                        rs->mR1PBeg = pbeg;
                        rs->mR1PEnd = pend;
                        rs->mR1SBeg = sbeg;
                        rs->mR1SEnd = send;
                        rs->mR1Seed = 1;
                    }
                }else{
                    rs->mR2SVID = svid;
                    rs->mR2MapQ = b->core.qual;
                    rs->mR2SRT = srst;
                    rs->mR2Hit = phit + shit;
                    if(saval){
                        rs->mR2PHit = phit;
                        rs->mR2SHit = shit;
                        rs->mR2PChr = pchr;
                        rs->mR2SChr = schr;
                        rs->mR2PBeg = pbeg;
                        rs->mR2PEnd = pend;
                        rs->mR2SBeg = sbeg;
                        rs->mR2Seed = 1;
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
                        iter->second->mR1PChr = pchr;
                        iter->second->mR1SChr = schr;
                        iter->second->mR1PBeg = pbeg;
                        iter->second->mR1PEnd = pend;
                        iter->second->mR1SBeg = sbeg;
                        iter->second->mR1SEnd = send;
                        iter->second->mR1Seed = 1;
                    }
                }else{
                    iter->second->mR2SVID = svid;
                    iter->second->mR2MapQ = b->core.qual;
                    iter->second->mR2SRT = srst;
                    iter->second->mR2Hit = phit + shit;
                    if(saval){
                        iter->second->mR2PHit = phit;
                        iter->second->mR2SHit = shit;
                        iter->second->mR2PChr = pchr;
                        iter->second->mR2SChr = schr;
                        iter->second->mR2PBeg = pbeg;
                        iter->second->mR2PEnd = pend;
                        iter->second->mR2SBeg = sbeg;
                        iter->second->mR2SEnd = send;
                        iter->second->mR2Seed = 1;
                    }
                }
            }
        }
    }
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
}

#endif
