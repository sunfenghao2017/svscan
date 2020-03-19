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
    int32_t mR1Hit = 0;
    int8_t mR1SRT = -1;
    int32_t mR2SVID = -1;
    uint8_t mR2MapQ = 0;
    int8_t mR2SRT = -1;
    int32_t mR2Hit = 0;
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
        int32_t lhit = 0, thit = 0;
        if(saval){
            int32_t sclen = 0;
            for(uint32_t i = 0; i < b->core.n_cigar; ++i){
                if(bam_cigar_op(bam_get_cigar(b)[i]) == BAM_CSOFT_CLIP){
                    sclen = bam_cigar_oplen(bam_get_cigar(b)[i]);
                    break;
                }
            }
            std::string qseq(b->core.l_qseq, '\0');
            for(int32_t i = 0; i < b->core.l_qseq; ++i){
                qseq[i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];
            }
            std::string lseq = qseq.substr(0, sclen);
            std::string tseq = qseq.substr(sclen);
            lhit = rf->validSRSeq(lseq);
            thit = rf->validSRSeq(tseq);
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
                }else{
                    rs->mR2SVID = svid;
                    rs->mR2MapQ = b->core.qual;
                    rs->mR2SRT = srst;
                    rs->mR2Hit = lhit + thit;
                }
                rssm[qname] = rs;
            }else{
                if(b->core.flag & BAM_FREAD1){
                    iter->second->mR1SVID = svid;
                    iter->second->mR1MapQ = b->core.qual;
                    iter->second->mR1SRT = srst;
                    iter->second->mR1Hit = lhit + thit;
                }else{
                    iter->second->mR2SVID = svid;
                    iter->second->mR2MapQ = b->core.qual;
                    iter->second->mR2SRT = srst;
                    iter->second->mR2Hit = lhit + thit;
                }
            }
        }
    }
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
}

#endif
