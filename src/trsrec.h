#ifndef TRSREC_H
#define TRSREC_H

#include <string>
#include <sstream>

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
    int32_t exon;        ///< exon number catenated
    int32_t insl;        ///< insertion length
    int32_t eoffset;     ///< offset to end of unit
    int32_t ioffset;     ///< offset to beg of unit
    bool pf5incc;        ///< 5' part of this unit catenated into fusion gene if true
    std::string mstat;   ///< matching status of this unit in catenation
    bool fullincc;       ///< this unit is all in cc
    std::string chr;     ///< used only in RNA sv

    /** TrsRec constructor */
    TrsRec(){
        name = "-";
        gene = "-";
        unit = "-";
        strand = "-";
        number = "-1";
        version = "-";
        primary = "-";
        exon = -1;
        pf5incc = true;
        fullincc = true;
    }

    /** TrsRec destructor */
    ~TrsRec(){}

    /** convert TrsRec to string
     * @return string representation of TrsRec
     */
    inline std::string toStr(){
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
        ret.append(",");
        ret.append(std::to_string(exon));
        return ret;
    }

    /** get name,uint,strand,number rep of TrsRec
     * @return string rep of TrsRec in name,strand,unit,number,exon format
     */
    inline std::string getTrs(){
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
        ret.append(",");
        ret.append(std::to_string(exon));
        return ret;
    }
    
    /** get name,uint,strand,number rep of TrsRec
     * @return string rep of TrsRec in name,strand,unit,number,exon format
     */
    inline std::string getTrsWithVer(){
        std::string ret;
        ret.append(name);
        ret.append(",");
        ret.append(strand);
        ret.append(",");
        ret.append(unit);
        ret.append(",");
        ret.append(number);
        ret.append(",");
        ret.append(std::to_string(exon));
        return ret;
    }

    /** test whether two transcript unit is near each other
     * @param other reference of another transcript record
     * @return true if they'are near each other
     */
    inline bool near(const TrsRec& other, float minRate = 0.8){
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

    /** parsing fusion part of transcript */
    inline void getCatPart(int32_t svt, bool svstart, bool rnamode = true){
        bool fwdstrand = (strand[0] == '+');
        if(rnamode) fwdstrand = true;
        int32_t svcat = svt;
        if(svcat >= 5) svcat -= 5;
        if(svstart){
            switch(svcat){
                case 0: case 2: // 5' part in cc
                    if(fwdstrand){
                        pf5incc = true;
                    }else{
                        pf5incc = false;
                    }
                    break;
                case 1: case 3: // 3' part in cc
                    if(fwdstrand){
                        pf5incc = false;
                    }else{
                        pf5incc = true;
                    }
                    break;
                default:
                    break;
            }
        }else{
            switch(svcat){
                case 0: case 3: // 5' part in cc
                    if(fwdstrand){
                        pf5incc = true;
                    }else{
                        pf5incc = false;
                    }
                    break;
                case 1: case 2: // 3' part in cc
                    if(fwdstrand){
                        pf5incc = false;
                    }else{
                        pf5incc = true;
                    }
                    break;
                default:
                    break;
            }
        }
    }

    /** get cigar */
    inline void getCigar(){
        std::stringstream ss;
        // e
        ss << "E" << exon << "(";
        // h or t
        if(pf5incc) ss << "H";
        else ss << "T";
        // m and i
        if(fullincc){
            ss << "AM" << insl << "I";
        }else{
            if(pf5incc){
                ss << ioffset;
            }else{
                ss << eoffset;
            }
            ss << "M" << insl << "I";
        }
        ss << ")";
        mstat = ss.str();
    }
};

/** class to store a list of TrsRec */
typedef std::vector<TrsRec> TrsRecList;

#endif
