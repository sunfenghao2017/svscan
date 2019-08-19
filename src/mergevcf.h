#ifndef MERGEVCF_H
#define MERGEVCF_H

#include <map>
#include <string>
#include <vector>
#include <cstdint>
#include <htslib/vcf.h>
#include "svutil.h"

/** interval representation of a SV */
struct Interval{
    int32_t mStart; ///< SV starting position
    int32_t mEnd;   ///< SV ending position
    int32_t mScore; ///< SV score

    /** Interval constructor
     * @param start SV starting position
     * @param end SV ending position
     * @param score SV score
     */
    Interval(int32_t start, int32_t end, int32_t score){
        mStart = start;
        mEnd = end;
        mScore = score;
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

/** type to store Interval of each contigs */
typedef std::vector<std::vector<Interval>> ContigIntervals;

/** type to store contig and id value pairs */
typedef std::map<std::string, int32_t> ContigMap;

/** class to store merge options */
struct VCFMerger{
    bool filterForPass; ///< keep PASS record if true
    bool filterForPrecise;///< keep PRECISE record if true
    int32_t chunksize;
    int32_t svcounter;
    int32_t bpoffset;
    int32_t minsize;
    int32_t maxsize;
    int32_t coverage;
    float recoverlap;
    float vaf;
    std::string outfile;
    std::vector<std::string> infiles;

    /** VCFMerger constructor */
    VCFMerger();

    /** VCFMerger destructor */
    ~VCFMerger();

    /** get unique tid for each contigs across samples
     * @param cmap map to store contig name ~ contig id value pairs
     */
    void getContigMap(ContigMap& cmap);

    /** get Interval representation of each SV of some type
     * @param ci reference of ContigIntervals to store SV on each contig
     * @param cmp contig name ~ contig id value pairs
     * @param svt sv type
     */   
    void getIntervals(ContigIntervals& ci, ContigMap& cmap, int32_t svt);

    /** merge Interval supporting same SV
     * @param ici referene of ContigIntervals to store input Intervals of SV
     * @param oci reference of ContigIntervals to store output Interval of SV
     * @param svt sv type
     */
    void mergeIntervals(ContigIntervals& ici, ContigIntervals& oci, int32_t svt);

};

#endif
