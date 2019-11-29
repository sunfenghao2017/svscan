#ifndef EXTRA_ANNO_H
#define EXTRA_ANNO_H

#include "cgranges.h"
#include "trsrec.h"
#include "util.h"
#include <fstream>
#include <string>
#include <vector>

class ExtraAnno{
    public:
        cgranges_t* extCrgs; ///< cgranges_t to store bed information
        TrsRecList trecs;    ///< gene names of each record

    ExtraAnno(){
        extCrgs = NULL;
    }

    ExtraAnno(const std::string& ebed){
        init(ebed);
    }

    ~ExtraAnno(){
        if(extCrgs) cr_destroy(extCrgs);
    }

    void init(const std::string& bedf);

    TrsRecList anno(const std::string& chr, int32_t beg, int32_t end);
};
#endif
