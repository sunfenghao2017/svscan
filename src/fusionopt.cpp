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
        bcf_get_info_string(hdr, rec, "CATT", &ct, &nct);
        int32_t recsvt = svutil::str2svt(ct, svt);
        int32_t svStart = rec->pos;
        int32_t svEnd = svStart + 1;
        if(bcf_get_info_int32(hdr, rec, "SVEND", &svend, &nsvend) > 0) svEnd = *svend;
        const char* chr1Name = bcf_hdr_id2name(hdr, rec->rid);
        bcf_get_info_string(hdr, rec, "CHREND", &chr2, &nchr2);
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
    if(!mExtraAnnoList.empty()) mExtraAnnotator.init(mExtraAnnoList);
    if(!mFsRptList.empty()) initFusionRptRange();
    if(!mGeneCrdList.empty()) initGeneCrdRange();
    if(!mHotPartnerList.empty()) initHotPartner();
    mInitialized = true;
}

bool FusionOptions::matchHotDirec(const std::string& hgene, const std::string& tgene){
    if(!mInitialized) init();
    auto hiter = m5Partners.find(hgene);
    auto titer = m3Partners.find(tgene);
    if(hiter != m5Partners.end() || titer != m3Partners.end()) return true;
    return false;
}

int FusionOptions::hasWhiteGene(const std::string& hgene, const std::string& tgene){
    if(!mInitialized) init();
    if(mWhiteGenes.empty()) return true;
    auto hiter = mWhiteGenes.find(hgene);
    auto titer = mWhiteGenes.find(tgene);
    int ret = 0;
    if(hiter != mWhiteGenes.end()) ret |= 1;
    if(titer != mWhiteGenes.end()) ret |= 2;
    return ret;
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
    auto rit = mBlackFusions.find(tgene);
    if(rit != mBlackFusions.end()){
        if(rit->second.find(hgene) != rit->second.end()){
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

void FusionOptions::initFusionRptRange(){
    if(mFsRptList.empty()) return;
    std::ifstream fr(mFsRptList);
    std::string line;
    std::vector<std::string> vstr;
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        std::string fsn = vstr[0] + "->" + vstr[2];
        if(vstr[2] == "*"){
            FusionRange sfsr;
            sfsr.mOneMatchOkay = true;
            mFusionRptMap[fsn] = sfsr;
        }else{
            if(vstr[1] == "*"){
                FusionRange tfsr;
                tfsr.mHgene = vstr[0];
                tfsr.mTgene = vstr[2];
                tfsr.mTwoMatchOkay = true;
                mFusionRptMap[fsn] = tfsr;
            }else{
                std::string uu = std::string(1, vstr[1][0]) + vstr[3][0];
                int32_t hu = std::atoi(vstr[1].substr(1).c_str());
                int32_t tu = std::atoi(vstr[3].substr(1).c_str());
                auto iter = mFusionRptMap.find(fsn);
                if(iter == mFusionRptMap.end()){
                    FusionRange fsr;
                    fsr.mHgene = vstr[0];
                    fsr.mTgene = vstr[2];
                    fsr.addp(hu, tu, uu);
                    mFusionRptMap[fsn] = fsr;
                }else{
                    iter->second.addp(hu, tu, uu);
                }
            }
        }
    }
    fr.close();
}

bool FusionOptions::inFsRptRange(std::string hgene, std::string tgene,  int32_t hu, int32_t tu, std::string uu){
    if(!mInitialized) init();
    std::string onemfs1 = hgene + "->*";
    if(mFusionRptMap.find(onemfs1) != mFusionRptMap.end()) return true;
    std::string onemfs2 = tgene + "->*";
    if(mFusionRptMap.find(onemfs2) != mFusionRptMap.end()) return true;
    std::string fgname = hgene + "->" + tgene;
    auto iter = mFusionRptMap.find(fgname);
    if(iter == mFusionRptMap.end()){
        return false;
    }else{
        if(iter->second.mTwoMatchOkay) return true;
        if(uu == "ee"){
            return iter->second.eegot(hu, tu);
        }else if(uu == "ei"){
            return iter->second.eigot(hu, tu);
        }else if(uu == "ii"){
            return iter->second.iigot(hu, tu);
        }else if(uu == "ie"){
            return iter->second.iegot(hu, tu);
        }
    }
    return false;
}

void FusionOptions::initGeneCrdRange(){
    if(mGeneCrdList.empty()) return;
    std::ifstream fr(mGeneCrdList);
    std::string line;
    std::vector<std::string> vstr;
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        GeneRange gr;
        gr.mGene = vstr[0];
        gr.mChr = vstr[1];
        gr.mStart = std::atoi(vstr[2].c_str());
        gr.mEnd = std::atoi(vstr[3].c_str());
        mGeneRangeVec.push_back(gr);
    }
    std::sort(mGeneRangeVec.begin(), mGeneRangeVec.end());
}
void FusionOptions::getInterGeneInfo(InterGeneInfo& igi, const std::string& chr, int32_t pos){
    if(mGeneRangeVec.empty()) return;
    if(!mInitialized) init();
    GeneRange gr;
    gr.mChr = chr;
    gr.mStart = pos;
    gr.mEnd = pos + 1;
    auto iter = std::lower_bound(mGeneRangeVec.begin(), mGeneRangeVec.end(), gr);
    if(iter == mGeneRangeVec.end() || iter->mChr != chr){ // reach end
        igi.dngene = "-";
        igi.dndist = 0;
    }else{
        igi.dngene = iter->mGene;
        igi.dndist = iter->mStart - pos;
    }
    if(iter == mGeneRangeVec.begin()){
        igi.upgene = "-";
        igi.updist = 0;
    }else{
        --iter;
        if(iter->mChr != chr){
            igi.upgene = "-";
            igi.updist = 0;
        }else{
            igi.upgene = iter->mGene;
            igi.updist = pos - iter->mEnd;
        }
    }
}

int32_t FusionOptions::geneNear(const std::string& g1, const std::string& chr1, int32_t pos1, const std::string& g2, const std::string& chr2){
    if(chr1 != chr2) return 0;
    if(g1 == g2) return 0;
    if(mGeneRangeVec.empty()) return 0;
    if(!mInitialized) init();
    GeneRange gr1;
    gr1.mChr = chr1;
    gr1.mStart = pos1;
    gr1.mEnd = pos1 + 1;
    auto iter = std::lower_bound(mGeneRangeVec.begin(), mGeneRangeVec.end(), gr1);
    while(iter != mGeneRangeVec.begin() && iter->mGene != g1){
        --iter;
    }
    if(iter->mGene != g1) return 0;
    gr1.mStart = iter->mStart;
    gr1.mEnd = iter->mEnd;
    if(iter != mGeneRangeVec.end()){
        ++iter;
        if(iter->mGene == g2){
            int32_t idst = iter->mStart - gr1.mEnd;
            if(idst == 0) return 1;
            else return idst; 
        }else{
            --iter;
        }
    }
    if(iter != mGeneRangeVec.begin()){
        --iter;
        if(iter->mGene == g2){
            int32_t idst = gr1.mStart - iter->mEnd;
            if(idst == 0) return 1;
            else return idst;
        }
    }
    return 0;
}

void FusionOptions::initHotPartner(){
    if(mHotPartnerList.empty()) return;
    std::ifstream fr(mHotPartnerList);
    std::string line;
    std::vector<std::string> vstr, hvs, ivs;
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        std::string fsgene = vstr[0] + "->" + vstr[1];
        HotPartner *hp = new HotPartner();
        util::split(vstr[2], hvs, ";");
        for(auto& e: hvs){
            util::split(e, ivs, ",");
            hp->hotpairs.push_back({std::atoi(ivs[0].c_str()), std::atoi(ivs[1].c_str())});
        }
        mHotPartnerMap[fsgene] = hp;
    }
}
