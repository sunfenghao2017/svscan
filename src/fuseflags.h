#ifndef FUSEFLAGS_H
#define FUSEFLAGS_H

#include "svutil.h"

struct FuseFlags{
    std::map<TFUSION_FLAG, std::string> flagmap;
    
    FuseFlags(){
        flagmap[FUSION_FALLGENE] = "partner are all genes";
        flagmap[FUSION_FNORMALCATDIRECT] = "catenation in normal pattern";
        flagmap[FUSION_FHOTGENE] = "partner in hotgene list";
        flagmap[FUSION_FCOMMONHOTDIRECT] = "hotgene partner in normal catenation pattern";
        flagmap[FUSION_FINDB] = "fusion pair in public db";
        flagmap[FUSION_FMIRROR] = "fusion has mirror event";
        flagmap[FUSION_FBLACKPAIR] = "fusion pair in blacklist";
        flagmap[FUSION_FFBG] = "fusion pair in background";
        flagmap[FUSION_FBLACKGENE] = "partner in blackgene list";
        flagmap[FUSION_FLOWCOMPLEX] = "partner has low complexity sequence";
        flagmap[FUSION_FPRIMARY] = "fusion should be reported as primary";
        flagmap[FUSION_FSUPPLEMENTARY] = "fusion should be reported as supplementary";
        flagmap[FUSION_FTOOSMALLSIZE] = "fusion size too small";
        flagmap[FUSION_FINSAMEGENE] = "fusion in same gene";
        flagmap[FUSION_FLOWAF] = "fusion af too low";
        flagmap[FUSION_FLOWSUPPORT] = "fusion support too small";
        flagmap[FUSION_FLOWDEPTH] = "fusion breakpoint depth too low";
        flagmap[FUSION_FPRECISE] = "fusion breakpoint is precise with concensus sequnce";
        flagmap[FUSION_FCALLFROMRNASEQ] = "fusion event is called from rna seq";
        flagmap[FUSION_FMIRRORINDB] = "fusion event's mirror event is in public db";
        flagmap[FUSION_FHTFLSWAPPED] = "fusion event's h/t gene swapped against bp1/2 gene";
        flagmap[FUSION_FREALNPASSED] = "fusion event has passed consensus sequence realignment";
        flagmap[FUSION_FINREPORTRNG] = "fusion event is in predefined report range";
    }

    ~FuseFlags(){}
    
    inline void showFuseFlags(){
        int32_t fwidth = std::numeric_limits<TFUSION_FLAG>::digits/4 + 4;
        for(auto& e: flagmap){
            std::cout << "  0x";
            std::cout <<  std::left << std::setw(fwidth);
            std::cout << std::hex << e.first << e.second << std::endl;
        }
    }
};

#endif
