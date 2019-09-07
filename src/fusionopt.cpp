#include "fusionopt.h"

FusionOptions::FusionOptions(){
    mWhiteFilter.mMinSupport = 3;
    mWhiteFilter.mMinVAF = 0;
    mUsualFilter.mMinSupport = 5;
    mUsualFilter.mMinVAF = 0.005;
};

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

void FusionOptions::parseFusionList(const std::string& fuseList, FusePairs& fusePairs){
    std::ifstream fr(fuseList);
    std::vector<std::string> vstr;
    std::string tmpstr;
    while(std::getline(fr, tmpstr)){
        util::split(tmpstr, vstr, "\t");
        fusePairs[vstr[0]].insert(vstr[3]);
    }
    fr.close();
}

void FusionOptions::getWhiteFusions(){
    parseFusionList(mWhiteList, mWhiteFusions);
}

void FusionOptions::getBlackFusions(){
    parseFusionList(mBlackList, mBlackFusions);
}

void FusionOptions::init(){
    mBgSVs.resize(9);
    if(!mBgBCF.empty()) getBgSVs();
    if(!mWhiteList.empty()) getWhiteFusions();
    if(!mBlackList.empty()) getBlackFusions();
    mInitialized = true;
}

bool FusionOptions::validFusion(const std::string& hgene, const std::string& tgene){
    if(!mInitialized) init();
    auto wit = mWhiteFusions.find(hgene);
    if(wit != mWhiteFusions.end()){
        if(wit->second.find(tgene) != wit->second.end()){
            return true;
        }
    }
    auto bit = mBlackFusions.find(hgene);
    if(bit != mBlackFusions.end()){
        if(bit->second.find(tgene) != bit->second.end()){
            return false;
        }
    }
    return true;
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
