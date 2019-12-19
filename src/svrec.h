#ifndef SVREC_H
#define SVREC_H

#include <string>
#include <cstdint>
#include <sstream>
#include <iostream>
#include "util.h"

struct SVRec{
    std::string svType;
    int32_t svSize;
    std::string bpMark;
    std::string bp1Chr;
    int32_t bp1Pos;
    std::string bp2Chr;
    int32_t bp2Pos;
    int32_t srCount;
    int32_t dpCount;
    int32_t srRescued;
    int32_t dpRescued;
    int32_t srRefCount;
    int32_t dpRefCount;
    int32_t molRescued;
    double af;
    int32_t insBp;
    std::string insSeq;
    std::string svSeq;
    int32_t seqBp;
    int32_t id;
    int32_t svInt;
    int32_t fsHits;
    std::string bp1Gene;
    std::string bp2Gene;
    std::string fuseGene;
    std::string fsMask;
    std::string trs1Name;
    int32_t trs1Pos;
    std::string trs2Name;
    int32_t trs2Pos;
    bool rnamode;

    SVRec(){
        rnamode = false;
    }

    ~SVRec(){}

    inline friend std::ostream& operator<<(std::ostream& os, const SVRec& svr){
        os << svr.svType << "\t" << svr.svSize << "\t" << svr.bpMark << "\t";
        os << svr.bp1Chr << "\t" << svr.bp1Pos << "\t" << svr.bp2Chr << "\t" << svr.bp2Pos << "\t";
        os << svr.srCount << "\t" << svr.dpCount << "\t";
        os << svr.srRescued << "\t" << svr.dpRescued << "\t" << svr.molRescued << "\t";
        os << svr.srRefCount << "\t" << svr.dpRefCount << "\t" << svr.af << "\t";
        os << svr.insBp << "\t" << svr.insSeq << "\t" << svr.svSeq << "\t" << svr.seqBp << "\t";
        os << svr.id << "\t" << svr.svInt << "\t";
        os << svr.bp1Gene << "\t" << svr.bp2Gene << "\t" << svr.fuseGene << "\t" << svr.fsMask << "\t" << svr.fsHits;
        if(svr.rnamode){
            os << "\t" << svr.trs1Name << "\t" << svr.trs1Pos << "\t" << svr.trs2Name << "\t" << svr.trs2Pos;
        }
        os << "\n";
        return os;
    }

    static std::string gethead(bool rnamode){
        std::stringstream oss;
        oss << "svType\tsvSize\tbpMark\t";//[0,2]
        oss << "bp1Chr\tbp1Pos\tbp2Chr\tbp2Pos\t"; // [3,6]
        oss << "srCount\tdpCount\t";// [7,8]
        oss << "srRescued\tdpRescued\tmolRescued\t"; // [9,11]
        oss << "srRefCount\tdpRefCount\tAF\t"; // [12,14]
        oss << "insBp\tinsSeq\tsvSeq\tseqBp\t";// [15,18]
        oss << "ID\tsvtInt\t"; // [19,20]
        oss << "bp1Gene\tbp2Gene\tfuseGene\tfsMask\tfsHits"; // [21,25]
        if(rnamode){
            oss << "\tts1Name\tts1Pos\tts2Name\tts2Pos\n"; //[26,29]
        }else{
            oss << "\n";
        }
        return oss.str();
    }

    static void line2rec(const std::string& line, SVRec& svr){
        std::vector<std::string> vstr;
        util::split(line, vstr, "\t");
        svr.svType = vstr[0];
        svr.svSize = std::atoi(vstr[1].c_str());
        svr.bpMark = vstr[2];
        svr.bp1Chr = vstr[3];
        svr.bp1Pos = std::atoi(vstr[4].c_str());
        svr.bp2Chr = vstr[5];
        svr.bp2Pos = std::atoi(vstr[6].c_str());
        svr.srCount = std::atoi(vstr[7].c_str());
        svr.dpCount = std::atoi(vstr[8].c_str());
        svr.srRescued = std::atoi(vstr[9].c_str());
        svr.dpRescued = std::atoi(vstr[10].c_str());
        svr.molRescued = std::atoi(vstr[11].c_str());
        svr.srRefCount = std::atoi(vstr[12].c_str());
        svr.dpRefCount = std::atoi(vstr[13].c_str());
        svr.af = std::atof(vstr[14].c_str());
        svr.insBp = std::atoi(vstr[15].c_str());
        svr.insSeq = vstr[16];
        svr.svSeq = vstr[17];
        svr.seqBp = std::atoi(vstr[18].c_str());
        svr.id = std::atoi(vstr[19].c_str());
        svr.svInt = std::atoi(vstr[20].c_str());
        svr.bp1Gene = vstr[21];
        svr.bp2Gene = vstr[22];
        svr.fuseGene = vstr[23];
        svr.fsMask = vstr[24];
        svr.fsHits = std::atoi(vstr[25].c_str());
        if(vstr.size() > 26){
            svr.trs1Name = vstr[26];
            svr.trs1Pos = std::atoi(vstr[27].c_str());
            svr.trs2Name = vstr[28];
            svr.trs2Pos = std::atoi(vstr[29].c_str());
        }
    }
};

#endif
