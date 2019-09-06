#ifndef INTERVAL_H
#define INTERVAL_H

#include <cstdint>
#include <iostream>

/** interval representation of a SV */
struct Interval{
    int32_t mStart; ///< SV starting position
    int32_t mEnd;   ///< SV ending position
    float mScore; ///< SV score

    /** Interval constructor
     * @param start SV starting position
     * @param end SV ending position
     */
    Interval(int32_t start, int32_t end){
        mStart = start;
        mEnd = end;
    }

    /** Interval constructor 
     * @param start SV starting position
     * @param end SV ending position
     * @param score SV score
     */
    Interval(int32_t start, int32_t end, float score){
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

#endif
