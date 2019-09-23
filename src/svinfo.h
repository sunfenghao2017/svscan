#ifndef SVINFO_H
#define SVINFO_H

#include <string>

/** Class to represent an SV info */
struct SVInfo{
    std::string mChr1; ///< large chr of an SV
    std::string mChr2; ///< lite chr of an SV
    int32_t mStart;    ///< SV starting position
    int32_t mEnd;      ///< SV ending position

    /** SVInfo constructor
     * @param chr1 large chr of an SV
     * @param chr2 lite chr of an SV
     * @param start SV starting position
     * @param end SV ending position
     */
    SVInfo(const std::string& chr1, const std::string& chr2, int32_t start, int32_t end){
        mChr1 = chr1;
        mChr2 = chr2;
        mStart = start;
        mEnd = end;
    }

    /** SVInfo destructor */
    ~SVInfo(){};
};

#endif
