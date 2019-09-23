#include "fusewlist.h"

void FuseWOpt::prepWlist(){
    // get genes from predefined fusion gene partner list
    util::loginfo("Beg parsing " + genelist + " to get gene list.");
    std::map<std::string, uint8_t> gset;
    std::ifstream fr(genelist);
    std::string tmpStr;
    std::vector<std::string> vstr;
    while(std::getline(fr, tmpStr)){
        util::split(tmpStr, vstr, "\t");
        uint8_t mask = 0;
        if(vstr[1].find("5") != std::string::npos) mask |= 2;
        if(vstr[1].find("3") != std::string::npos) mask |= 1;
        gset.insert({vstr[0], mask});
    }
    fr.close();
    util::loginfo("End parsing " + genelist + ", got " + std::to_string(gset.size()) + " genes ");
    // parse fusedb to get whitelist as subset
    std::ofstream fw(whitelist);
    fr.open(fusedb.c_str());
    std::string hgene, tgene;
    while(std::getline(fr, tmpStr)){
        util::split(tmpStr, vstr, "\t");
        hgene = vstr[0];
        tgene = vstr[1];
        bool hotMark[2] = {false, false};
        auto hiter = gset.find(hgene);
        if(hiter != gset.end()){
            hotMark[0] = true;
            hiter->second |= 4;
        }
        auto titer = gset.find(tgene);
        if(titer != gset.end()){
            hotMark[1] = true;
            titer->second |= 4;
        }
        if(hotMark[0] || hotMark[1]){
            std::string hotStr;
            if(hotMark[0]) hotStr.append("Y\t");
            else hotStr.append("N\t");
            if(hotMark[1]) hotStr.append("Y\t");
            else hotStr.append("N\t");
            if(hiter != gset.end()){
                if(hiter->second & 2) hotStr.append("\tY");
                else hotStr.append("\tN");
            }else hotStr.append("\t-");
            if(titer != gset.end()){
                if(titer->second & 1) hotStr.append("\tY");
                else hotStr.append("\tN");
            }else hotStr.append("\t-");
            fw << hgene << "\t" << tgene << "\t" << hotStr << "\n";
        }
    }
    // print warning information about hot gene not found in db
    for(auto& e: gset){
        if(!(e.second & 4)){
            std::cout << "Hot gene missed in out file: " << e.first << std::endl;
        }
    }
    fr.close();
    fw.close();
}
