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
        fw << FusionRecord::gethead(false, fromrna);
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
            // get exon-exon groups
            std::map<std::string, uint32_t> eecnt;
            for(uint32_t i = 0; i < iter->second.size(); ++i){
                std::string exstr = std::to_string(iter->second[i].exon1) + "," + std::to_string(iter->second[i].exon2);
                auto eit = eecnt.find(exstr);
                if(eit == eecnt.end()){
                    eecnt[exstr] = i;
                }else{
                    if(iter->second[i].fusionmols > iter->second[eit->second].fusionmols && iter->second[i].fuserate > iter->second[eit->second].fuserate && iter->second[i].srcount > 0){
                        eit->second = i;
                    }
                }
                iter->second[i].report = false;
            }
            for(auto& e: eecnt){
                iter->second[e.second].report = true;
            }
            // check whether primary reported
            int pridx = -1;
            for(uint32_t i = 0; i < iter->second.size(); ++i){
                if((iter->second[i].fsmask & FUSION_FPRIMARY) && iter->second[i].report){
                    pridx = i;
                }
            }
            if(pridx >= 0){// drop all others
                for(uint32_t i = 0; i < iter->second.size(); ++i){
                    if(iter->second[i].report && !(iter->second[i].fsmask & FUSION_FPRIMARY)){
                        iter->second[i].report = false;
                    }
                }
            }
        }
        // output
        for(auto iter= fm.begin(); iter != fm.end(); ++iter){
            for(uint32_t i = 0; i < iter->second.size(); ++i){
                if(iter->second[i].report){
                    // remove exon info
                    iter->second[i].transcript1 = iter->second[i].transcript1.substr(0, iter->second[i].transcript1.find_last_of(","));
                    iter->second[i].transcript2 = iter->second[i].transcript2.substr(0, iter->second[i].transcript2.find_last_of(","));
                    iter->second[i].outurl = false;
                    fw << iter->second[i];
                }
            }
        }
        // close file
        fw.close();
    }
};

#endif
