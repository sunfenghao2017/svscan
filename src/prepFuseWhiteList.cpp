#include <fstream>
#include <iostream>
#include <util.h>
#include <set>
#include <map>

int main(int argc, char** argv){
    if(argc < 4){
        std::cout << argv[0] << " <gene.tsv> <fusedb.tsv> <white.tsv>" << std::endl;
        return 0;
    }
    std::string gsetf = argv[1];
    std::string gfsdb = argv[2];
    std::string outf = argv[3];
    // get genes from predefined fusion gene partner list
    util::loginfo("Beg parsing " + gsetf + " to get gene list.");
    std::map<std::string, bool> gset;
    std::ifstream fr(gsetf);
    std::string tmpStr;
    while(std::getline(fr, tmpStr)){
        gset.insert({tmpStr, false});
    }
    fr.close();
    util::loginfo("End parsing " + gsetf + ", got " + std::to_string(gset.size()) + " genes ");
    // parse fusedb to get whitelist as subset
    std::ofstream fw(outf);
    fr.open(gfsdb.c_str());
    std::vector<std::string> vstr;
    std::string hgene, tgene;
    while(std::getline(fr, tmpStr)){
        util::split(tmpStr, vstr, "\t");
        hgene = vstr[0];
        tgene = vstr[1];
        bool hotMark[2] = {false, false};
        auto hiter = gset.find(hgene);
        if(hiter != gset.end()){
            hotMark[0] = true;
            hiter->second = true;
        }
        auto titer = gset.find(tgene);
        if(titer != gset.end()){
            hotMark[1] = true;
            titer->second = true;
        }
        if(hotMark[0] || hotMark[1]){
            std::string hotStr;
            if(hotMark[0]) hotStr.append("Y\t");
            else hotStr.append("N\t");
            if(hotMark[1]) hotStr.append("Y\n");
            else hotStr.append("N\n");
            fw << hgene << "\t" << tgene << "\t" << hotStr;
        }
    }
    // print warning information about hot gene not found in db
    for(auto& e: gset){
        if(!e.second){
            std::cout << "Hot gene missed in out file: " << e.first << std::endl;
        }
    }
    fr.close();
    fw.close();
}
