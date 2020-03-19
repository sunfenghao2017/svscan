#ifndef REALNFILTER_H
#define REALNFILTER_H

#include "onlinebwa.h"
#include <bamutil.h>

struct BpPair{
    int32_t tid1;
    int32_t tid2;
    int32_t pos1;
    int32_t pos2;
    bool swapped = false;

    BpPair(){}

    ~BpPair(){}

    void adjustpt();

    bool agree(const BpPair& other);
};

class RealnFilter{
    public:
        OnlineBWA* mBWA;
        bam_hdr_t* mHeader;

    RealnFilter(){
        mBWA = NULL;
        mHeader = NULL;
    }

    RealnFilter(const std::string& ref){
        init(ref);
    }

    ~RealnFilter(){
        if(mBWA) delete mBWA;
        if(mHeader) bam_hdr_destroy(mHeader);
    }

    void init(const std::string& ref){
        mBWA = new OnlineBWA();
        mBWA->loadIndex(ref);
        mHeader =mBWA->getBamHeader();
    }

    /** realignment test of concensus sequence
     * return integer more than 4 if perfect secondary pairs > 2 
     * return 0 if nice perfect match or just one primary sc
     * return -1 if whole seq match continuous
     * return -2 if two primary sc match
     * return -3 if bp not match
     */
    int32_t validCCSeq(const std::string& seq, const std::string& chr1, int32_t& pos1, const std::string& chr2, int32_t& pos2, int32_t fseq, int32_t inslen);

    /** realignment test of split read sequence */
    int32_t validSRSeq(const std::string& seq);
};

#endif
