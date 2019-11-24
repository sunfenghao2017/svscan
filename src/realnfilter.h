#ifndef REALNFILTER_H
#define REALNFILTER_H

#include "onlinebwa.h"
#include "bamutil.h"

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

    bool validCCSeq(const std::string& seq, const std::string& chr1, int32_t& pos1, const std::string& chr2, int32_t& pos2);
};

#endif
