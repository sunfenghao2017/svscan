#include "svfilter.h"

SVFilter::SVFilter(){
    softEnv = new Software();
    softEnv->cmp += "version: " + softEnv->ver + "\n";
    softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
    bgSVs.resize(9);
}

SVFilter::~SVFilter(){
    if(softEnv) delete softEnv;
}

void SVFilter::update(int argc, char** argv){
    // update software environment records
    softEnv->cwd = util::cwd();
    for(int i = 0; i < argc; ++i){
        softEnv->cmd.append(argv[i]);
        softEnv->cmd.append(" ");
    }
}

void SVFilter::getBgSVs(){
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
        bgSVs[recsvt].push_back({chr1Name, chr2, svStart, svEnd});
    }
    if(svt) free(svt);
    if(ct) free(ct);
    if(chr2) free(chr2);
    if(svend) free(svend);
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
}

void SVFilter::init(){
    if(!initialized){
        getBgSVs();
        initialized = true;
    }
}

bool SVFilter::validSV(int32_t svt, const std::string& chr1, const std::string& chr2, int32_t start, int32_t end){
    if(!initialized) init();
    for(auto& e: bgSVs[svt]){
        if(e.mChr1 == chr1 && e.mChr2 == chr2){
            if(std::abs(e.mStart - start) < maxBpOffset && std::abs(e.mEnd - end) < maxBpOffset){
                return false;
            }
        }
    }
    return true;
}

void SVFilter::filter(){
    init();
    std::ifstream fr(infile);
    std::string tmpstr;
    std::vector<std::string> vstr;
    std::getline(fr, tmpstr);
    std::ofstream fw(outfile);
    fw << tmpstr << "\n";
    while(std::getline(fr, tmpstr)){
        util::split(tmpstr, vstr, "\t");
        int32_t svt = std::atoi(vstr[29].c_str());
        int32_t start = std::stoi(vstr[11].c_str());
        int32_t end = std::atoi(vstr[14].c_str());
        std::string chr1 = vstr[10];
        std::string chr2 = vstr[13];
        int32_t sr = std::atoi(vstr[18].c_str());
        int32_t dp = std::atoi(vstr[19].c_str());
        float af = std::atof(vstr[22].c_str());
        if(sr < minsr) continue;
        if(dp < mindp) continue;
        if(af < minaf) continue;
        if(!validSV(svt, chr1, chr2, start, end)) continue;
        fw << tmpstr << "\n";
    }
    fr.close();
    fw.close();
}
