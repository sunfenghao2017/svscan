#include "region.h"

Region::Region(const std::string& regfile, const std::string& bamfile){
    samFile* fp = sam_open(bamfile.c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(fp);
    std::ifstream fr(regfile);
    std::vector<std::string> vstr;
    std::string tmpStr;
    std::vector<std::vector<std::pair<int32_t, int32_t>>> oriReg(hdr->n_targets);
    while(std::getline(fr, tmpStr)){
        util::split(tmpStr, vstr, "\t");
        int32_t tid = bam_name2id(hdr, vstr[0].c_str());
        oriReg[tid].push_back(std::make_pair(std::atoi(vstr[1].c_str()), std::atoi(vstr[2].c_str())));
    }
    for(int32_t i = 0; i < hdr->n_targets; ++i){
        std::vector<std::pair<int32_t, int32_t>> ctgReg = mergeAndSortRegions(oriReg[i]);
    }
    for(int32_t i = 0; i < hdr->n_targets; ++i){
        for(auto& e: oriReg[i]){
            mRegs[i].insert(e);
        }
    }
}

std::vector<std::pair<int32_t, int32_t>> Region::mergeAndSortRegions(std::vector<std::pair<int32_t, int32_t>>& regs){
    std::vector<std::pair<int32_t, int32_t>> ret;
    std::stack<std::pair<int32_t, int32_t>> s;
    std::sort(regs.begin(),
              regs.end(),
              [](std::pair<int32_t, int32_t>& p1, std::pair<int32_t, int32_t>& p2){
                 return p1.first < p2.first || (p1.first == p2.first && p1.second < p2.second);});
    s.push(regs[0]);
    for(uint32_t i = 1; i < regs.size(); ++i){
        std::pair<int32_t, int32_t> lastReg = s.top();
        std::pair<int32_t, int32_t> curReg = regs[i];
        if(curReg.first >= lastReg.first && curReg.first <= lastReg.second){
            std::pair<int32_t, int32_t> newReg = std::make_pair(lastReg.first, std::max(lastReg.second, curReg.second));
            s.pop();
            s.push(newReg);
        }else{
            s.push(curReg);
        }
    }
    while(!s.empty()){
        ret.push_back(s.top());
        s.pop();
    }
    std::sort(ret.begin(),
          ret.end(),
          [](std::pair<int32_t, int32_t>& p1, std::pair<int32_t, int32_t>& p2){
             return p1.first < p2.first || (p1.first == p2.first && p1.second < p2.second);});
    return ret;
}

std::set<std::pair<int32_t, int32_t>>::iterator Region::getFirstOverlap(int32_t tid, const std::pair<int32_t, int32_t>& p){
    auto iter= std::upper_bound(mRegs[tid].begin(), mRegs[tid].end(), p);
    auto upIter = iter;
    if(iter != mRegs[tid].begin()){
        auto lowIter = --iter;
        bool lowGet = false;
        while(p.first <= lowIter->second){
            if(iter == mRegs[tid].begin()){
                return iter;
            }
            --lowIter;
            lowGet = true;
        }
        if(lowGet){
            return ++lowIter;
        }
    }
    bool upGet = false;
    while(upIter != mRegs[tid].end() && p.second >= upIter->first){
        ++upIter;
        upGet = true;
    }
    if(upGet){
        return --upIter;
    }
    return mRegs[tid].end();
}
