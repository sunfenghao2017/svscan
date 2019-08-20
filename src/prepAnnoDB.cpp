#include <map>
#include <util.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>

int main(int argc, char** argv){
    if(argc < 3){
        printf("%s <refGene.Acc.tsv.gz> <refGene.Anno.gz>\n", argv[0]);
        return 0;
    }
    BGZF* ifp = bgzf_open(argv[1], "rb");
    BGZF* ofp = bgzf_open(argv[2], "wb");
    kstring_t str = {0, 0, 0};
    bgzf_getline(ifp, '\n', &str);
    std::vector<std::string> vstr;
    std::vector<std::string> istr, estr;
    std::vector<int32_t> iint, eint;
    std::stringstream oss;
    // chr start end strand feature count trsname tgenename trsversion
    while(bgzf_getline(ifp, '\n', &str) > 0){
        oss.str(""); 
        oss.clear();
        util::split(str.s, vstr, "\t");
        std::string trs = vstr[1];
        std::string chr = vstr[2];
        std::string gene = vstr[12];
        std::string version = vstr[16];
        std::string strand = vstr[3];
        int32_t trsStart = std::atoi(vstr[4].c_str());
        int32_t trsEnd = std::atoi(vstr[5].c_str()) - 1;
        int32_t cdsStart = std::atoi(vstr[6].c_str());
        int32_t cdsEnd = std::atoi(vstr[7].c_str()) - 1;
        // utr range
        int32_t utr5start = trsStart;
        int32_t utr5end = cdsStart - 1;
        int32_t utr3start = cdsEnd + 1;
        int32_t utr3end = trsEnd;
        if(strand[0] == '-'){
            utr5start = cdsEnd + 1;
            utr5end = trsEnd;
            utr3start = trsStart;
            utr3end = cdsStart - 1;
        }
        if(utr5start < utr5end){
            oss << chr << "\t" << utr5start << "\t" << utr5end << "\t" << strand << "\tutr5\t0\t" << trs << "\t" << gene << "\t" << version << "\n";
        }
        if(utr3start < utr3end){
            oss << chr << "\t" << utr3start << "\t" << utr3end << "\t" << strand << "\tutr3\t0\t" << trs << "\t" << gene << "\t" << version << "\n";
        }
        // exon range
        util::split(vstr[9], istr, ",");
        util::split(vstr[10], estr, ",");
        util::strvec2intvec(istr, iint);
        util::strvec2intvec(estr, eint);
        for(uint32_t i = 0; i < iint.size(); ++i){
            oss << chr << "\t" << iint[i] << "\t" << (eint[i] - 1) << "\t" << strand << "\texon\t";
            if(strand[0] == '+') oss << i + 1;
            else oss << (iint.size() - i);
            oss << "\t" << trs << "\t" << gene << "\t" << version << "\n";
        }
        // intron range
        std::vector<int32_t> intronbeg;
        std::vector<int32_t> intronend;
        for(uint32_t i = 0; i < iint.size() - 1; ++i){
            if(eint[i] < iint[i + 1] - 1){
                intronbeg.push_back(eint[i]);
                intronend.push_back(iint[i + 1] - 1);
            }
        }
        for(uint32_t i = 0; i < intronbeg.size(); ++i){
            oss << chr << "\t" << intronbeg[i] << "\t" << intronend[i] << "\t" << strand << "\tintron\t";
            if(strand[0] == '+') oss << i + 1;
            else oss << (intronbeg.size() - i);
            oss << "\t" << trs << "\t" << gene << "\t" << version << "\n";
        }
        // output str
        std::string recs = oss.str();
        assert(bgzf_write(ofp, recs.c_str(), recs.size()) >= 0);
    }
    bgzf_close(ofp);
    bgzf_close(ifp);
}
