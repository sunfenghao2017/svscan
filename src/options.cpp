#include "options.h"

Options::Options(){
    madCutoff = 9;
    nthread = 8;
    bcfOut = "out.bcf";
    filterOpt = new SVFilter();
    softEnv = new Software();
    msaOpt = new MSAOpt();
    passOpt = new PassOptions();
    softEnv->cmp += "version: " + softEnv->version + "\n";
    softEnv->cmp += "updated: " + std::string(__TIME__);
    libInfo = NULL;
    contigNum = 0;
}

Options::~Options(){
    if(filterOpt) delete filterOpt;
    if(softEnv) delete softEnv;
    if(libInfo) delete libInfo;
}

void Options::validate(){
    // check genome indexed
    if(!util::exists(genome + ".fai")){
        util::errorExit("Genome must be indexed by `samtools faidx`");
    }
}

void Options::update(int argc, char** argv){
    // update software environment records
    softEnv->cwd = util::cwd();
    for(int i = 0; i < argc; ++i){
        softEnv->cmd.append(argv[i]);
        softEnv->cmd.append(" ");
    }
    // get library information
    libInfo = getLibInfo(bamfile);
    // update SV types to discover
    std::vector<std::string> svt = {"INV", "DEL", "DUP", "INS", "BND"};
    std::string allSVT;
    if(svtypes.size() == 0){
        for(int i = 0; i < 9; ++i) SVTSet.insert(i);
        allSVT.append("INV DEL DUP INS BND");
    }else{
        for(auto& e: svtypes){
            allSVT.append(svt[e] + " ");
            if(e == 0){// Inversion
                SVTSet.insert(0);
                SVTSet.insert(1);
            }else if(e == 4){// translocation
                SVTSet.insert({5, 6, 7, 8});
            }else{// Deletion, Duplication, Insertion
                SVTSet.insert(e);
            }
        }
    }
    // update contigNum
    contigNum = libInfo->mContigNum;
    // show SV types to discover
    util::loginfo("SV types to discover: " + allSVT);
}

LibraryInfo* Options::getLibInfo(const std::string& bam){
    LibraryInfo* libInfo = new LibraryInfo();
    samFile* fp = sam_open(bam.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    std::vector<int32_t> vecISize, readLen;
    const uint16_t BAM_SKIP_RECORD_MASK = (BAM_FREAD2 | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP);
    while(sam_read1(fp, h, b) >= 0){
        if(b->core.flag & BAM_SKIP_RECORD_MASK) continue;
        readLen.push_back(b->core.l_qseq);
        if(b->core.flag & BAM_FPAIRED && b->core.tid == b->core.mtid){
            vecISize.push_back(std::abs(b->core.isize));
        }
    }
    libInfo->mReadLen = statutil::median(readLen);
    libInfo->mMedian = statutil::median(vecISize);
    libInfo->mMad = statutil::mad(vecISize, libInfo->mMedian);
    libInfo->mMaxNormalISize = libInfo->mMedian + (5 * libInfo->mMad);
    libInfo->mMinNormalISize = std::max(0, libInfo->mMedian - (5 * libInfo->mMad));
    libInfo->mMinISizeCutoff = std::max(0, libInfo->mMedian - (madCutoff * libInfo->mMad));
    libInfo->mMaxISizeCutoff = std::max(libInfo->mMedian + (madCutoff * libInfo->mMad), 2 * libInfo->mReadLen);
    libInfo->mMaxISizeCutoff = std::max(libInfo->mMaxISizeCutoff, 500);
    libInfo->mContigNum = h->n_targets;
    libInfo->mVarisize = std::max(libInfo->mReadLen, libInfo->mMaxNormalISize);
    sam_close(fp);
    bam_destroy1(b);
    bam_hdr_destroy(h);
    return libInfo;
}

void Options::getValidRegion(){
    samFile* fp = sam_open(bamfile.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    validRegions.resize(h->n_targets);
    // Parse invalid region if exists
    RegionList invalidRegions;
    invalidRegions.resize(h->n_targets);
    if(!excludeReg.empty()){
        std::ifstream fr(excludeReg);
        std::vector<std::string> vstr;
        std::string tmpStr;
        while(std::getline(fr, tmpStr)){
            util::split(tmpStr, vstr, "\t");
            int32_t tid = bam_name2id(h, vstr[0].c_str());
            invalidRegions[tid].insert({std::atoi(vstr[1].c_str()), std::atoi(vstr[2].c_str())});
        }
    }
    // Create valid regions
    for(int32_t i = 0; i < h->n_targets; ++i){
        int32_t istart = 0;
        for(auto iter = invalidRegions[i].begin(); iter != invalidRegions[i].end(); ++iter){
            if(istart < iter->first) validRegions[i].insert({istart, iter->first - 1});
            istart = iter->second + 1;
        }
        if(istart < (int32_t)h->target_len[i]) validRegions[i].insert({istart, (int32_t)h->target_len[i] - 1});
    }
    // Clean-up
    sam_close(fp);
    bam_hdr_destroy(h);
}
