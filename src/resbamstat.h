#ifndef RES_BAM_STAT_H
#define RES_BAM_STAT_H

#include <htslib/sam.h>
#include <cstdint>
#include <string>
#include <map>

/** class to store a read supporting status */
struct ReadSupport{
    int32_t mR1SVID = -1;
    uint8_t mR1MapQ = 0;
    int8_t mR1SRT = -1;
    int32_t mR2SVID = -1;
    uint8_t mR2MapQ = 0;
    int8_t mR2SRT = -1;
};

/** type to store read supporting statistics */
typedef std::map<std::string, ReadSupport*> ReadSupportStatMap;

inline void getReadSupportStatus(const std::string& bam, ReadSupportStatMap& rssm){
    samFile* fp = sam_open(bam.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    while(sam_read1(fp, h, b) >= 0){
        uint8_t* idata = bam_aux_get(b, "ZF");
        uint8_t* sdata = bam_aux_get(b, "ST");
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
                }else{
                    rs->mR2SVID = svid;
                    rs->mR2MapQ = b->core.qual;
                    rs->mR2SRT = srst;
                }
                rssm[qname] = rs;
            }else{
                if(b->core.flag & BAM_FREAD1){
                    iter->second->mR1SVID = svid;
                    iter->second->mR1MapQ = b->core.qual;
                    iter->second->mR1SRT = srst;
                }else{
                    iter->second->mR2SVID = svid;
                    iter->second->mR2MapQ = b->core.qual;
                    iter->second->mR2SRT = srst;
                }
            }
        }
    }
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
}

#endif
