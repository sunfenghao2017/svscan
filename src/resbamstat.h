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
    int32_t mR1SHit = 0;
    std::string mR1SChr;
    int32_t mR1SBeg = -1;
    int32_t mR1SEnd = -1;
    
    int32_t mR2SVID = -1;
    uint8_t mR2MapQ = 0;
    int8_t mR2SRT = -1;
    int32_t mR2Hit = 0;
    int32_t mR2PHit = 0;
    std::string mR2PChr;
    int32_t mR2PBeg = -1;
    int32_t mR2PEnd = -1;
    int32_t mR2SHit = 0;
    std::string mR2SChr;
    int32_t mR2SBeg = -1;
    int32_t mR2SEnd = -1;
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
        int32_t lhit = 0, lbeg = -1, lend = -1, thit = 0, tbeg = -1, tend = -1;
        std::string lchr, tchr;
        if(saval){
            int32_t sclen = 0;
            for(uint32_t i = 0; i < b->core.n_cigar; ++i){
                if(bam_cigar_op(bam_get_cigar(b)[i]) == BAM_CSOFT_CLIP){
                    sclen = bam_cigar_oplen(bam_get_cigar(b)[i]);
                    break;
                }
            }
            lchr = h->target_name[b->core.tid];
            lbeg = b->core.pos;
            lend = bam_endpos(b);
            std::string qseq(b->core.l_qseq, '\0');
            for(int32_t i = 0; i < b->core.l_qseq; ++i){
                qseq[i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];
            }
            std::string lseq = qseq.substr(0, sclen);
            std::string tseq = qseq.substr(sclen);
            lhit = rf->validSRSeq(lseq);
            thit = rf->validSRSeq(tseq);
            std::string ss = bam_aux2Z(saval);
            std::vector<std::string> vsas;
            util::split(ss, vsas, ";");
            std::vector<std::string> vstr;
            util::split(vsas[0], vstr, ",");
            tchr = vstr[0];
            tbeg = std::atoi(vstr[1].c_str());
            tend = tbeg + sclen;
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
                    rs->mR1Hit = lhit + thit;
                    if(saval){
                        rs->mR1PHit = lhit;
                        rs->mR1SHit = thit;
                        rs->mR1PChr = lchr;
                        rs->mR1SChr = tchr;
                        rs->mR1PBeg = lbeg;
                        rs->mR1PEnd = lend;
                        rs->mR1SBeg = tbeg;
                        rs->mR1SEnd = tend;
                    }
                }else{
                    rs->mR2SVID = svid;
                    rs->mR2MapQ = b->core.qual;
                    rs->mR2SRT = srst;
                    rs->mR2Hit = lhit + thit;
                    if(saval){
                        rs->mR2PHit = lhit;
                        rs->mR2SHit = thit;
                        rs->mR2PChr = lchr;
                        rs->mR2SChr = tchr;
                        rs->mR2PBeg = lbeg;
                        rs->mR2PEnd = lend;
                        rs->mR2SBeg = tbeg;
                    }
                }
                rssm[qname] = rs;
            }else{
                if(b->core.flag & BAM_FREAD1){
                    iter->second->mR1SVID = svid;
                    iter->second->mR1MapQ = b->core.qual;
                    iter->second->mR1SRT = srst;
                    iter->second->mR1Hit = lhit + thit;
                    if(saval){
                        iter->second->mR1PHit = lhit;
                        iter->second->mR1SHit = thit;
                        iter->second->mR1PChr = lchr;
                        iter->second->mR1SChr = tchr;
                        iter->second->mR1PBeg = lbeg;
                        iter->second->mR1PEnd = lend;
                        iter->second->mR1SBeg = tbeg;
                        iter->second->mR1SEnd = tend;
                    }
                }else{
                    iter->second->mR2SVID = svid;
                    iter->second->mR2MapQ = b->core.qual;
                    iter->second->mR2SRT = srst;
                    iter->second->mR2Hit = lhit + thit;
                    if(saval){
                        iter->second->mR2PHit = lhit;
                        iter->second->mR2SHit = thit;
                        iter->second->mR2PChr = lchr;
                        iter->second->mR2SChr = tchr;
                        iter->second->mR2PBeg = lbeg;
                        iter->second->mR2PEnd = lend;
                        iter->second->mR2SBeg = tbeg;
                        iter->second->mR2SEnd = tend;
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
