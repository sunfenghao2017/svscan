#ifndef SVREC_H
#define SVREC_H

#include <string>
#include <cstdint>
#include <sstream>
#include <iostream>

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
        os << svr.srRescued << "\t" << svr.dpRescued << "\t";
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
        oss << "srRescued\tdpRescued\t"; // [9,10]
        oss << "srRefCount\tdpRefCount\tAF\t"; // [11,13]
        oss << "insBp\tinsSeq\tsvSeq\tseqBp\t";// [14,17]
        oss << "ID\tsvtInt\t"; // [18,19]
        oss << "bp1Gene\tbp2Gene\tfuseGene\tfsMask\tfsHits"; // [20,24]
        if(rnamode){
            oss << "\tts1Name\tts1Pos\tts2Name\tts2Pos\n"; //[25,28]
        }else{
            oss << "\n";
        }
        return oss.str();
    }
};

#endif
