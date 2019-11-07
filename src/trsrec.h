#ifndef TRSREC_H
#define TRSREC_H

#include <string>

/** class to store an transcript record */
struct TrsRec{
    std::string name;    ///< transcript name
    std::string gene;    ///< gene name
    std::string unit;    ///< uint name
    std::string strand;  ///< strand
    std::string number;  ///< count
    std::string version; ///< transcript version
    std::string primary; ///< Y if it's a canonical transcript
    bool drop;           ///< drop from report if true
    int32_t pos;         ///< position on genome
    int32_t eoffset;     ///< offset to end of unit  
    int32_t ioffset;     ///< offset to beg of unit
    std::string chr;     ///< used only in RNA sv

    /** TrsRec constructor */
    TrsRec(){}

    /** TrsRec destructor */
    ~TrsRec(){}

    /** convert TrsRec to string
     * @return string representation of TrsRec
     */
    std::string toStr(){
        std::string ret;
        ret.append(gene);
        ret.append(",");
        ret.append(name);
        ret.append(".");
        ret.append(version);
        ret.append(",");
        ret.append(unit);
        ret.append(",");
        ret.append(strand);
        ret.append(",");
        ret.append(number);
        return ret;
    }

    /** get name,uint,strand,number rep of TrsRec
     * @return string rep of TrsRec in name,strand,unit,number format
     */
    std::string getTrs(){
        std::string ret;
        ret.append(name);
        ret.append(".");
        ret.append(version);
        ret.append(",");
        ret.append(strand);
        ret.append(",");
        ret.append(unit);
        ret.append(",");
        ret.append(number);
        return ret;
    }

    /** test whether two transcript unit is near each other
     * @param other reference of another transcript record
     * @return true if they'are near each other
     */
    bool near(const TrsRec& other, float minRate = 0.8){
        if(name != other.name) return false; // not same transcript
        if(unit == other.unit && number == other.number){
            if(unit == "intron"){
                return true; // all in same intron
            }else if(unit == "exon"){
                int32_t exsize = eoffset + ioffset;
                return std::abs(ioffset - other.ioffset) / (double)exsize < minRate;
            }else{
                return true; // utr region
            }
        }
        if(unit != other.unit){
            if(std::abs(std::atoi(number.c_str()) - std::atoi(other.number.c_str())) <= 1){// nearby exon/intron
                return true;
            }
        }
        return false;
    }
};

/** class to store a list of TrsRec */
typedef std::vector<TrsRec> TrsRecList;

#endif
