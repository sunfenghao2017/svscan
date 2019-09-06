#include "filter.h"

Filter::Filter(){
    softEnv = new Software();
    softEnv->cmp += "version: " + softEnv->version + "\n";
    softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
}

Filter::~Filter(){
    if(softEnv){
        delete softEnv;
    }
}

void Filter::update(int argc, char** argv){
    // update software environment records
    softEnv->cwd = util::cwd();
    for(int i = 0; i < argc; ++i){
        softEnv->cmd.append(argv[i]);
        softEnv->cmd.append(" ");
    }
}

void Filter::getBkgSVIntervals(SVIntervals& si){
    si.resize(9);
    htsFile* fp = bcf_open(bkgbcf.c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    bcf1_t* rec = bcf_init();
    char* svt = NULL;
    int32_t nsvt = 0;
    char* ct = NULL;
    int32_t nct = 0;
    int32_t* svend = NULL;
    int32_t nsvend = 0;
    while(bcf_read(fp, hdr, rec) >= 0){
        bcf_unpack(rec, BCF_UN_INFO);
        bcf_get_info_string(hdr, rec, "SVTYPE", &svt, &nsvt);
        bcf_get_info_string(hdr, rec, "CT", &ct, &nct);
        int32_t recsvt = svutil::str2svt(ct, svt);
        int32_t svStart = rec->pos;
        int32_t svEnd = svStart + 1;
        if(bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend) > 0) svEnd = *svend;
        const char* chrName = bcf_hdr_id2name(hdr, rec->rid);
        si[recsvt][chrName].push_back({svStart, svEnd});
    }
    if(svt) free(svt);
    if(ct) free(ct);
    if(svend) free(svend);
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
}

void Filter::parseFusePairs(FusePairs& fusePairs, const std::string& fuseFile){
    std::ifstream fr(fuseFile);
    std::vector<std::string> vstr;
    std::string tmpstr;
    while(std::getline(fr, tmpstr)){
        util::split(tmpstr, vstr, "\t");
        fusePairs[vstr[0]].push_back(vstr[3]);
    }
    fr.close();
}

int main(int argc, char** argv){

}
