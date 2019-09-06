#ifndef FILTER_H
#define FILTER_H

#include <map>
#include <set>
#include <string>
#include <htslib/vcf.h>
#include "svutil.h"
#include "util.h"

/** class to store software environment */
struct Software{
    std::string version = "0.0.0"; ///< software version
    std::string cmd;               ///< software execution command
    std::string cwd;               ///< software execution directory
    std::string cmp;               ///< software update time

    /** Software constructor */
    Software(){}

    /** Software destructor */
    ~Software(){}
};

/** interval representation of a SV */
struct Interval{
    int32_t mStart; ///< SV starting position
    int32_t mEnd;   ///< SV ending position

    /** Interval constructor
     * @param start SV starting position
     * @param end SV ending position
     */
    Interval(int32_t start, int32_t end){
        mStart = start;
        mEnd = end;
    }

    /** Interval destructor */
    ~Interval(){}

    /** operator to compare two Interval object
     * @param other reference of another Interval object
     * @return true if this Interval < other
     */
    inline bool operator<(const Interval& other) const {
        return mStart < other.mStart || (mStart == other.mStart && mEnd < other.mEnd);
    }

    /** operator to output interval to ostream
     * @param os reference of ostream
     * @param i reference of Interval
     * @return reference of os
     */
    inline friend std::ostream& operator<<(std::ostream& os, const Interval& i){
        os << "==================Interval==================\n";
        os << "Starting Pos: " << i.mStart << "\n";
        os << "Ending Pos: " << i.mEnd << "\n";
        os << "============================================\n";
        return os;
    }

    /** get overlap ratio of two Interval object
     * @param other reference of another Interval object
     * @return overlap ratio of these two Interval objects
     */
    inline double overlap(const Interval& other) const {
        if(mEnd < other.mStart || other.mEnd < mStart) return 0;
        double lenA = mEnd - mStart;
        if(lenA <= 0) return 0;
        double lenB = other.mEnd - other.mStart;
        if(lenB <= 0) return 0;
        double overlapLen = std::min(mEnd, other.mEnd) - std::max(mStart, other.mStart);
        if(overlapLen <= 0) return 0;
        return overlapLen/std::max(lenA, lenB);
    }
};

/** type to store one kind of structural variant Intervals on each contigs */
typedef std::map<std::string, std::vector<Interval>> ContigIntervals;

/** type to store each kind of structural variant Intervals */
typedef std::vector<ContigIntervals> SVIntervals;

/** type to stroe fusion gene partner information hgene ~ tgene list */
typedef std::map<std::string, std::vector<std::string>> FusePairs;

/** class to do filter of a bcf file containing structural variants */
struct Filter{
    bool passFilter = false;    ///< keep PASS variant if true
    bool preciseFilter = false; ///< keep PRECISE variant if true
    int32_t minsize = 0;        ///< min variant size to keep
    int32_t maxsize = (1>>20);  ///< max variant size to keep
    int32_t mincov = 0;         ///< min coverage of variant to keep
    float minvaf = 0;           ///< min VAF of variant to keep
    std::string whitelist;      ///< white list of structural variants
    std::string blacklist;      ///< black list of structural variants
    std::string bkgbcf;         ///< structural variants in control sample(bcf formats)
    std::string outfile;        ///< output file
    std::string infile;         ///< input file
    Software* softEnv;         ///< pointer to Software object
    
    /** default constructor of Filter */
    Filter();
    
    /** default destructor of Filter */
    ~Filter();
   
    /** update options after command arguments parsed
     * @param argc number of arguments feed to main
     * @param argv array of arguments feed to main
     */
    void update(int argc, char** argv);

    /** get Interval representation of all SVs in bkgbcf
     * @param si reference of SVIntervals to store all SVs on all contigs
     */
    void getBkgSVIntervals(SVIntervals& si);

    /** parse fusion gene partner infomation table
     * @param fusePairs reference of FusePairs
     * @param fuseFile fusion partner database file
     */
    void parseFusePairs(FusePairs& fusePairs, const std::string& fuseFile);

    /** do filter job */
    void filter();
};

#endif
