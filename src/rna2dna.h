#ifndef RNA2DNA_H
#define RNA2DNA_H

#include <string>
#include <cstdint>
#include <util.h>

/** transcript unit to dna coordinate record */
struct Rna2DnaUnit{
    std::string tname;    ///< transcript name
    int32_t tbeg;         ///< transcript beginning position
    int32_t tend;         ///< transcript ending position
    std::string uname;    ///< unit name
    int32_t ucount;       ///< unit count
    std::string gname;    ///< gene name
    std::string gchr;     ///< chromosome
    int32_t gbeg;         ///< genome beginning position
    int32_t gend;         ///< genome ending position
    char gstrand;         ///< genome strand
    std::string tversion; ///< transcript version

    /** default constructor */
    Rna2DnaUnit(){}

    /** construct from a record of rna annodb */
    Rna2DnaUnit(const std::string& rec){
        std::vector<std::string> vstr;
        util::split(rec, vstr, "\t");
        tname = vstr[0];
        tbeg = std::atoi(vstr[1].c_str());
        tend = std::atoi(vstr[2].c_str());
        uname = vstr[3];
        ucount = std::atoi(vstr[4].c_str());
        gname = vstr[5];
        gchr = vstr[6];
        gbeg = std::atoi(vstr[7].c_str());
        gend = std::atoi(vstr[8].c_str());
        gstrand = vstr[9][0];
        tversion = vstr[10];
    }

    /** default destructor */
    ~Rna2DnaUnit(){}

    /** operator to compare two Rna2DnaUnit */
    inline bool operator<(const Rna2DnaUnit& other) const {
        return ucount < other.ucount;
    }

    /** whether an pos of this transcript is included by this unit */
    bool incpos(const int32_t& bp){
        return bp <= tend && bp >= tbeg;
    } 

    /** test whether this unit containing the bp is the proper catenated exon
     * @param bp breakpoint position
     * @param svt sv type
     * @param start if true bp is bp1
     * @return 0 if this unit is proper, positive if next one is proper, negative if previous if proper, return value is inslen value
     */
    std::pair<int32_t, int32_t> catthis(const int32_t& bp, int32_t svt, bool start, int moff = 10){
        if(svt >= 5) svt -= 5;
        switch(svt){
            case 0: // 5t5 catenation
                if(bp - tbeg < moff){
                    return {tbeg - bp, -1};
                }
                break;
            case 1: // 3t3 catenation
                if(tend - bp < moff){
                    return {tend - bp, +1};
                }
                break;
            case 2: // 5t3 catenation
                if(start){
                    if(bp - tbeg < moff){
                        return {tbeg - bp, -1};
                    }
                }else{
                    if(tend - bp < moff){
                        return {tend - bp, +1};
                    }
                }
                break;
            case 3: // 3t5 catenation
                if(start){
                    if(tend - bp < moff){
                        return {tend - bp, +1};
                    }
                }else{
                    if(bp - tbeg < moff){
                        return {tbeg - bp, -1};
                    }
                }
                break;
            default:
                break;
        }
        return {0, 0};
    }
};

#endif
