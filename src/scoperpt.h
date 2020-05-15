#ifndef SCRPT_H
#define SCRPT_H

#include "fuserec.h"

struct ScopeRptOpt{
    std::string infs;
    std::string outfs;
    bool fromrna = false;

    void out(){
        // write out head
        std::ofstream fw(outfs);
        fw << FusionRecord::gethead(fromrna);
        // collect fs
        std::ifstream fr(infs);
        std::string line;
        std::getline(fr, line);
        FusionRecord f;
        std::map<std::string, std::vector<FusionRecord>> fm;
        FusionRecSchema s;
        while(std::getline(fr, line)){
            f.line2rec(line, s);
            fm[f.fusegene].push_back(f);
        }
        fr.close();
        // choose to output
        for(auto iter= fm.begin(); iter != fm.end(); ++iter){
            // get highest af one
            double maxaf = iter->second[0].fuserate;
            int maxfr = iter->second[0].fusionmols;
            uint32_t maxai = 0;
            for(uint32_t i = 0; i < iter->second.size(); ++i){
                if(iter->second[i].fuserate > maxaf && iter->second[i].fusionmols > maxfr){
                    maxaf = iter->second[i].fuserate;
                    maxfr = iter->second[i].fusionmols;
                    maxai = i;
                }
            }
            for(uint32_t i = 0; i < iter->second.size(); ++i){
                if(i != maxai) iter->second[i].report = false;
                else iter->second[i].report = true;
            }
        }
        // output
        for(auto iter= fm.begin(); iter != fm.end(); ++iter){
            for(uint32_t i = 0; i < iter->second.size(); ++i){
                if(iter->second[i].report){
                    // remove intron info
                    iter->second[i].transcript1 = iter->second[i].transcript1.substr(0, iter->second[i].transcript1.find_first_of(","));
                    iter->second[i].transcript2 = iter->second[i].transcript2.substr(0, iter->second[i].transcript2.find_first_of(","));
                    fw << iter->second[i];
                }
            }
        }
        // close file
        fw.close();
    }
};

#endif
