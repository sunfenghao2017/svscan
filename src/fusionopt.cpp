#include "fusionopt.h"

void FusionOptions::getBgSVs(){
    mBgSVs.resize(9);
    htsFile* fp = bcf_open(mBgBCF.c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    bcf1_t* rec = bcf_init();
    char* svt = NULL;
    int32_t nsvt = 0;
    char* ct = NULL;
    int32_t nct = 0;
    int32_t* svend = NULL;
    int32_t nsvend = 0;
    char* chr2 = NULL;
    int32_t nchr2 = 0;
    while(bcf_read(fp, hdr, rec) >= 0){
        bcf_unpack(rec, BCF_UN_INFO);
        bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
        bcf_get_info_string(hdr, rec, "CT", &ct, &nct);
        int32_t recsvt = svutil::str2svt(ct, svt);
        int32_t svStart = rec->pos;
        int32_t svEnd = svStart + 1;
        if(bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) svEnd = *svend;
        const char* chr1Name = bcf_hdr_id2name(hdr, rec->rid);
        bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2);
        mBgSVs[recsvt].push_back({chr1Name, chr2, svStart, svEnd});
    }
    if(svt) free(svt);
    if(ct) free(ct);
    if(chr2) free(chr2);
    if(svend) free(svend);
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
}

void FusionOptions::parseWhiteList(){
    std::ifstream fr(mWhiteList);
    std::vector<std::string> vstr;
    std::string tmpstr;
    while(std::getline(fr, tmpstr)){
        util::split(tmpstr, vstr, "\t");
        mWhiteFusions[vstr[0]].insert(vstr[1]);
        if(vstr[2] == "Y") mWhiteGenes.insert(vstr[0]);
        if(vstr[3] == "Y") mWhiteGenes.insert(vstr[1]);
        if(vstr[4] == "Y") m5Partners.insert(vstr[0]);
        if(vstr[5] == "Y") m3Partners.insert(vstr[1]);
    }
    fr.close();
}

void FusionOptions::parseBlackList(){
    std::ifstream fr(mBlackList);
    std::vector<std::string> vstr;
    std::string tmpstr;
    while(std::getline(fr, tmpstr)){
        util::split(tmpstr, vstr, "\t");
        if(vstr[1] != "*") mBlackFusions[vstr[0]].insert(vstr[1]);
        else mBlackGenes.insert(vstr[0]);
    }
    fr.close();
}

void FusionOptions::parseSameGeneEventList(){
    std::ifstream fr(mSameGeneSVList);
    std::vector<std::string> vstr;
    std::vector<std::string> sve;
    std::vector<std::string> svt;
    std::string tmpstr;
    while(std::getline(fr, tmpstr)){
        util::split(tmpstr, vstr, "\t");
        util::split(vstr[1], sve, ",");
        util::split(vstr[2], svt, ",");
        DetectRange dr;
        dr.mGene = vstr[0];
        if(sve[0] != "."){
            for(auto& e: sve) dr.mExonList.insert(std::atoi(e.c_str()));
        }
        for(auto& t: svt) dr.mSVT.insert(std::atoi(t.c_str()));
        mDetectRngMap[dr.mGene] = dr;
    }
    fr.close();
}
void FusionOptions::init(){
    mBgSVs.resize(9);
    if(!mBgBCF.empty()) getBgSVs();
    if(!mWhiteList.empty()){
        parseWhiteList();
    }
    if(!mBlackList.empty()) parseBlackList();
    if(!mSameGeneSVList.empty()) parseSameGeneEventList();
    mInitialized = true;
}

bool FusionOptions::matchHotDirec(const std::string& hgene, const std::string& tgene){
    if(!mInitialized) init();
    auto hiter = m5Partners.find(hgene);
    auto titer = m3Partners.find(tgene);
    if(hiter != m5Partners.end() || titer != m3Partners.end()) return true;
    return false;
}

bool FusionOptions::hasWhiteGene(const std::string& hgene, const std::string& tgene){
    if(!mInitialized) init();
    auto hiter = mWhiteGenes.find(hgene);
    auto titer = mWhiteGenes.find(tgene);
    if(hiter != mWhiteGenes.end() || titer != mWhiteGenes.end()) return true;
    return false;
}

bool FusionOptions::hasBlackGene(const std::string& hgene, const std::string& tgene){
    if(!mInitialized) init();
    auto hiter = mBlackGenes.find(hgene);
    auto titer = mBlackGenes.find(tgene);
    if(hiter != mBlackGenes.end() || titer != mBlackGenes.end()) return true;
    return false;
}

bool FusionOptions::inWhiteList(const std::string& hgene, const std::string& tgene){
    if(!mInitialized) init();
    auto wit = mWhiteFusions.find(hgene);
    if(wit != mWhiteFusions.end()){
        if(wit->second.find(tgene) != wit->second.end()){
            return true;
        }
    }
    return false;
}

bool FusionOptions::inBlackList(const std::string& hgene, const std::string& tgene){
    auto bit = mBlackFusions.find(hgene);
    if(bit != mBlackFusions.end()){
        if(bit->second.find(tgene) != bit->second.end()){
            return true;
        }
    }
    return false;
}

bool FusionOptions::inSameSVRngMap(const std::string& gene, const std::vector<int32_t>& exon, int32_t svt){
    auto giter = mDetectRngMap.find(gene);
    if(giter == mDetectRngMap.end()) return false;
    auto titer = giter->second.mSVT.find(svt);
    if(titer == giter->second.mSVT.end()) return false;
    if(!giter->second.mExonList.empty()){
        for(auto& e: exon){
            if(giter->second.mExonList.find(e) != giter->second.mExonList.end()) return true;
        }
        return false;
    }else{
        return true;
    }
}

bool FusionOptions::validSV(int32_t svt, const std::string& chr1, const std::string& chr2, int32_t start, int32_t end){
    if(!mInitialized) init();
    for(auto& e: mBgSVs[svt]){
        if(e.mChr1 == chr1 && e.mChr2 == chr2){
            if(std::abs(e.mStart - start) < mMaxBpOffset && std::abs(e.mEnd - end) < mMaxBpOffset){
                return false;
            }
        }
    }
    return true;
}
