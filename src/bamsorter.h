#ifndef BAMSORTER_H
#define BAMSORTER_H

#include <htslib/sam.h>

struct BamComp{
    inline bool operator()(const bam1_t* b1, const bam1_t* b2) const {
        if(b1->core.tid >= 0) {        // b1 is mapped
            if(b2->core.tid<0 )
                return true;
            else if(b2->core.tid >  b1->core.tid)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos >  b1->core.pos)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid >  b1->core.mtid)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid == b1->core.mtid && b2->core.mpos >  b1->core.mpos)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid == b1->core.mtid && b2->core.mpos == b1->core.mpos) {
                if(b2->core.isize > b1->core.isize)
                    return true;
                else
                    return (long)b2->data > (long)b1->data;
            } else
                return false;
        } else {         // b1 is unmapped
            if(b2->core.tid<0) { // both are unmapped
                return (long)b2->data > (long)b1->data;
            }
            else
                return false;
        }
    }
};

#endif
