#ifndef PROD_RPT_H
#define PROD_RPT_H

#include "fuserec.h"

struct ProdRptOpt{
    std::string infs;
    std::string outfp;
    std::string outfm;
    std::string hlist;
    bool fromrna = false;
    std::set<std::string> hotgene;

    void get_hot(){
        std::ifstream fr(hlist);
        std::vector<std::string> vstr;
        std::string tmpstr;
        while(std::getline(fr, tmpstr)){
            util::split(tmpstr, vstr, "\t");
            if(vstr[2] == "Y") hotgene.insert(vstr[0]);
            if(vstr[3] == "Y") hotgene.insert(vstr[1]);
        }
        fr.close();
    }

    void out(){
        // write out head
        std::ofstream fwp(outfp);
        fwp << FusionRecord::gethead(true, fromrna);
        std::ofstream fwm(outfm);
        fwm << FusionRecord::gethead(true, fromrna);
        // collect fs
        std::ifstream fr(infs);
        std::string line;
        std::getline(fr, line);
        FusionRecord f;
        std::map<std::string, std::vector<FusionRecord>> fm; //hotgene->fs
        FusionRecSchema s;
        while(std::getline(fr, line)){
            f.line2rec(line, s);
            f.jctpos1 -= 1;
            f.jctpos2 -= 1;
            auto iter = hotgene.find(f.gene1);
            if(iter != hotgene.end()) fm[f.gene1].push_back(f);
            else{
                iter = hotgene.find(f.gene2);
                if(iter != hotgene.end()) fm[f.gene2].push_back(f);
            }
        }
        fr.close();
        // choose to output
        for(auto iter= fm.begin(); iter != fm.end(); ++iter){
            // test YY
            std::set<uint32_t> yyf;
            for(uint32_t i = 0; i < iter->second.size(); ++i){
                if(iter->second[i].indb == "Y" && iter->second[i].status.find_first_of("CD") == std::string::npos) yyf.insert(i);
            }
            if(yyf.empty()){// no yy, output all
               for(uint32_t i = 0; i < iter->second.size(); ++i){
                  iter->second[i].report = true;
               }
            }else{// only output yy
                for(uint32_t i = 0; i < iter->second.size(); ++i){
                    if(yyf.find(i) == yyf.end()) iter->second[i].report = false;
                    else iter->second[i].report = true;
                }
            }
        }
        // output
        for(auto iter= fm.begin(); iter != fm.end(); ++iter){
            for(uint32_t i = 0; i < iter->second.size(); ++i){
                // remove exon info
                iter->second[i].transcript1 = iter->second[i].transcript1.substr(0, iter->second[i].transcript1.find_last_of(","));
                iter->second[i].transcript2 = iter->second[i].transcript2.substr(0, iter->second[i].transcript2.find_last_of(","));
                iter->second[i].outurl = true;
                // fix YM to YY
                if(iter->second[i].indb == "Y" && iter->second[i].status == "M"){
                    iter->second[i].status = "Y";
                }
                iter->second[i].outrec(fwm);
                if(iter->second[i].report) iter->second[i].outrec(fwp);
            }
        }
        // close file
        fwp.close();
        fwm.close();
    }
};

#endif
