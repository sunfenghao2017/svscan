#include "svdnadb.h"

void SVDNADBOpt::prepDB(){
    // parse gene2main transcript map
    std::map<std::string, std::string> g2cmap;
    util::makeMapPairFromFileByLine(gene2cnc, g2cmap);
    BGZF* ifp = bgzf_open(refSeqDB.c_str(), "rb");
    BGZF* ofp = bgzf_open(svAnnoDB.c_str(), "wb");
    std::ofstream ofc(bedCDS); // cds bed region
    std::ofstream ofu(bedUnits); // unit bed region
    std::ofstream ofr(ref2gene); // ref2gene tsv
    kstring_t str = {0, 0, 0};
    bgzf_getline(ifp, '\n', &str);
    std::vector<std::string> vstr;
    std::vector<std::string> istr, estr;
    std::vector<int32_t> iint, eint;
    std::stringstream aoss; // annotation db record buffer
    std::stringstream coss; // cds bed region record buffer
    std::stringstream uoss; // unit bed region buffer
    std::stringstream ross; // ref2gene buffer
    // chr start end strand feature count trsname tgenename trsversion
    while(bgzf_getline(ifp, '\n', &str) > 0){
        aoss.str(""); coss.str(""); uoss.str(""); ross.str("");
        aoss.clear(); coss.clear(); uoss.clear(); ross.clear();
        util::split(str.s, vstr, "\t");
        std::string trs = vstr[1];
        std::string chr = vstr[2];
        std::string gene = vstr[12];
        std::string version = vstr[16];
        std::string strand = vstr[3];
        ross << trs << "\t" << gene << "\n";
        int32_t trsStart = std::atoi(vstr[4].c_str());
        int32_t trsEnd = std::atoi(vstr[5].c_str()) - 1;
        int32_t cdsStart = std::atoi(vstr[6].c_str());
        int32_t cdsEnd = std::atoi(vstr[7].c_str()) - 1;
        std::string msMK = "Y"; // canonical transcript marker
        auto iter = g2cmap.find(gene);
        if(iter != g2cmap.end()){ // has a canonical transcript predefiend
            if(trs != iter->second) msMK = "N";
        }
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
            aoss << chr << "\t" << utr5start << "\t" << utr5end << "\t" << strand << "\tutr5\t0\t" << trs << "\t" << gene << "\t" << version << "\t" << msMK << "\n";
            uoss << chr << "\t" << utr5start << "\t" << utr5end + 1 << "\t" << trs << "\tutr5\t";
            uoss << gene << "\t" << strand << "\t" << version << "\t" << msMK << "\n";
        }
        if(utr3start < utr3end){
            aoss << chr << "\t" << utr3start << "\t" << utr3end << "\t" << strand << "\tutr3\t0\t" << trs << "\t" << gene << "\t" << version << "\t" << msMK << "\n";
            uoss << chr << "\t" << utr3start << "\t" << utr3end + 1 << "\t" << trs << "\tutr3\t";
            uoss << gene << "\t" << strand << "\t" << version << "\t" << msMK << "\n";
        }
        // exon range
        util::split(util::rstrip(vstr[9], ","), istr, ",");
        util::split(util::rstrip(vstr[10], ","), estr, ",");
        util::strvec2intvec(istr, iint);
        util::strvec2intvec(estr, eint);
        for(uint32_t i = 0; i < iint.size(); ++i){
            aoss << chr << "\t" << iint[i] << "\t" << (eint[i] - 1) << "\t" << strand << "\texon\t";
            uoss << chr << "\t" << iint[i] << "\t" << eint[i] << "\t" << trs << "\texon";
            if(cdsStart < cdsEnd){
                if(!(iint[i] > cdsEnd) && !(eint[i] - 1 < cdsStart)){
                    coss << chr << "\t" << std::max(iint[i], cdsStart) << "\t" << std::min(eint[i], cdsEnd + 1) << "\t" << trs << "\t" << "exon\t";
                    if(strand[0] == '+') coss << i + 1 << "\t" << gene << "\t" << strand << "\t" << version << "\t" << msMK << "\n";
                    else coss << (iint.size() - i) << "\t" << gene << "\t" << strand << "\t" << version << "\t" << msMK << "\n";
                }
            }
            if(strand[0] == '+'){
                aoss << i + 1;
                uoss << i + 1;
            }
            else{
                aoss << (iint.size() - i);
                uoss << (iint.size() - i);
            }
            aoss << "\t" << trs << "\t" << gene << "\t" << version << "\t" << msMK << "\n";
            uoss << "\t" << gene << "\t" << strand << "\t" << version << "\t" << msMK << "\n";
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
            aoss << chr << "\t" << intronbeg[i] << "\t" << intronend[i] << "\t" << strand << "\tintron\t";
            uoss << chr << "\t" << intronbeg[i] << "\t" << intronend[i] + 1 << "\t" << trs << "\tintron";
            if(strand[0] == '+'){
                aoss << i + 1;
                uoss << i + 1;
            }else{
                aoss << (intronbeg.size() - i);
                uoss << (intronbeg.size() - i);
            }
            aoss << "\t" << trs << "\t" << gene << "\t" << version << "\t" << msMK << "\n";
            uoss << "\t" << gene << "\t" << strand << "\t" << version << "\t" << msMK << "\n";
        }
        // output str
        std::string recs = aoss.str();
        ofc << coss.str();
        ofu << uoss.str();
        ofr << ross.str();
        assert(bgzf_write(ofp, recs.c_str(), recs.size()) >= 0);
    }
    bgzf_close(ofp);
    bgzf_close(ifp);
    ofc.close();
    ofu.close();
    ofr.close();
}
