#ifndef ALNDESCRIPTOR_H
#define ALNDESCRIPTOR_H

#include <cstdint>
#include <iostream>

/** alignment result descriptor of SR consensus sequence against SV ref sequence */
struct AlignDescriptor{
    int32_t mCSStart = 0;        ///< gap starting position on consensus sequence(0-based), next alignment is gap
    int32_t mCSEnd = 0;          ///< gap ending position on consensus sequence(1-based), this base is first after gap
    int32_t mRefStart = 0;       ///< gap starting position on SV ref sequence(0-based), next alignment is gap
    int32_t mRefEnd = 0;         ///< gap ending position on SV ref sequence(1-based), this base is first after gap
    int32_t mHomLeft = 0;        ///< longest homology length between left part of gapped seq and the other seq ended at gap end position(reversely) 
    int32_t mHomRight = 0;       ///< longest homology length between right part of gapped seq and the other seq started at gap start position
    double mPercID = 0;          ///< alignment identity percents of this alignment skipping leading tailing and this continuous gap
    int32_t mMinGapNeeded = 15;  ///< minimal gap length needed for an valid alignment
    int32_t mMaxCoordOffset = 5; ///< maximum coordinate offset allowed for gapped seq
    double mMinPercID = 0.90;    ///< minimal alignment identity percents needed for an valid alignment

    /** Constructor */
    AlignDescriptor() = default;
    
    /** Destructor */
    ~AlignDescriptor() = default;

    /** Operator to output an AlignDescriptor to ostream
     * @param os reference of ostream object
     * @param ad reference of AlignDescriptor object
     * @return reference of ostream
     */
    inline friend std::ostream& operator<<(std::ostream& os, const AlignDescriptor& ad){
        os << "===========================================================================================================================\n";
        os << "Gap starting position on consensus sequence(0-based): " << ad.mCSStart << "\n";
        os << "Gap ending position on consensus sequence(1-based): " << ad.mCSEnd << "\n";
        os << "Gap starting position on constructed SV reference sequence(0-based): " << ad.mRefStart << "\n";
        os << "Gap ending position on constructed SV reference sequence(1-based): " << ad.mRefEnd << "\n";
        os << "Out of gap alignment identity percentage: " << ad.mPercID << "\n";
        os << "Predefined minimal inner gap length needed: " << ad.mMinGapNeeded << "\n";
        os << "longest homology length between left part of gapped seq and the other seq ended at gap end position(reversely): " << ad.mHomLeft << "\n";
        os << "ongest homology length between right part of gapped seq and the other seq started at gap start position: " << ad.mHomRight << "\n";
        os << "Predefined maximum starting/ending gap coordinates offset on gapped sequence: " << ad.mMaxCoordOffset << "\n";
        os << "Predefined minimal alignment identity percentage of gap alignment skipping leading, tailing and this continuous gap: " << ad.mMinPercID << "\n";
        os << "===========================================================================================================================\n";
        return os;
    }

    /** check whether this gapped alignment is valid
     * @return true if valid
     */
    inline bool validGapAlignment(int32_t svt){
        if(svt == 4) return (mRefEnd - mRefStart < mMaxCoordOffset) && (mCSEnd - mCSStart > mMinGapNeeded);
        return (mRefEnd - mRefStart > mMinGapNeeded) && (mCSEnd - mCSStart < mMaxCoordOffset);
    }
    
    /** check whether this gapped alignment have an valid identity ratio
     * @return true if valid
     */
    inline bool validFlankQual(){
        return mPercID >= mMinPercID;
    }
};

#endif
