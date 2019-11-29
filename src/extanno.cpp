#include "extanno.h"

void ExtraAnno::init(const std::string& bedf){
    extCrgs = cr_init();
    std::ifstream fr(bedf);
    std::string line;
    std::vector<std::string> vstr;
    int32_t index = 0;
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        cr_add(extCrgs, vstr[0].c_str(), std::atoi(vstr[1].c_str()), std::atoi(vstr[2].c_str()), index);
        TrsRec tr;
        tr.chr = vstr[0];
        tr.gene = vstr[3];
        tr.name = vstr[4];
        trecs.push_back(tr);
        ++index;
    }
    cr_index(extCrgs);
}

TrsRecList ExtraAnno::anno(const std::string& chr, int32_t beg, int32_t end){
    TrsRecList trl;
    int64_t n, *b = 0, max_b = 0;
    n = cr_overlap(extCrgs, chr.c_str(), beg, end, &b, &max_b);
    for(int64_t i = 0; i < n; ++i){
        trl.push_back(trecs[cr_label(extCrgs, b[i])]);
    }
    free(b);
    return trl;
}
