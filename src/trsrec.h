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

};

/** class to store a list of TrsRec */
typedef std::vector<TrsRec> TrsRecList;

#endif
