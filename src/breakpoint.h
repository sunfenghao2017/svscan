#ifndef BREAKPOINT_H
#define BREAKPOINT_H

#include <cstdint>
#include <string>
#include <htslib/sam.h>
#include "svrecord.h"

/** strategy to construct pseudo reference
 * F stands for forward strand of reference 
 * R stands for reverse complement of reference
 * bp stands for breakpoint
 * svt: 0 ------------F-----------bp------------R-----------
 * svt: 1 ------------R-----------bp------------F-----------
 * svt: 2 ------------F-----------bp------------F-----------
 * svt: 3 ------------F-----------bp------------F-----------
 * svt: 4 ------------F-----------bp------------F-----------
 * svt: 5 ------------F-----------bp------------R-----------
 * svt: 6 ------------R-----------bp------------F-----------
 * svt: 7 ------------F-----------bp------------F-----------
 * svt: 8 ------------F-----------bp------------F-----------
 */

/** class to store and analysis breakpoint of an SV */
class BreakPoint{
    public:
        int32_t mSVStartBeg = 0; ///< SV leftmost starting position
        int32_t mSVStartEnd = 0; ///< SV rightmost starting position
        int32_t mSVEndBeg = 0;   ///< SV leftmost ending position
        int32_t mSVEndEnd = 0;   ///< SV rightmost ending position
        int32_t mSVStart = 0;    ///< SV starting position
        int32_t mSVEnd = 0;      ///< SV ending position
        int32_t mPESupport = 0;  ///< SV pair-end support number
        int32_t mSRSupport = 0;  ///< SV split-read support number
        int32_t mSVT = 0;        ///< SV type
        int32_t mChr1 = 0;       ///< (larger) chr of SV
        int32_t mChr2 = 0;       ///< (little) chr of SV
        int32_t mChr1Len = 0;    ///< (larger) chr length of SV
        int32_t mChr2Len = 0;    ///< (little) chr length of SV
        int32_t mBoundary = 0;   ///< boundary of start and end position of SV(in SR, it's the consensus seq of SRs)

    public:

        /** BreakPoint default constructor */
        BreakPoint(){}

        /** Construct a BreakPoint from an existing SV 
         * @param sv pointer of SVRecord 
         * @param hdr bam header
         */
        BreakPoint(const SVRecord* sv, const bam_hdr_t* hdr){
            mSVStartBeg = sv->mSVStart;
            mSVStartEnd = sv->mSVStart;
            mSVEndBeg = sv->mSVEnd;
            mSVEndEnd = sv->mSVEnd;
            mSVStart = sv->mSVStart;
            mSVEnd = sv->mSVEnd;
            mPESupport = sv->mPESupport;
            mSRSupport = sv->mSRSupport;
            mSVT = sv->mSVT;
            mChr1 = sv->mChr1;
            mChr2 = sv->mChr2;
            mChr1Len = hdr->target_len[mChr1];
            mChr2Len = hdr->target_len[mChr2];
            mBoundary = sv->mConsensus.size();
            if(hdr) init(mChr1Len, mChr2Len);
        }

        /** Construct a BreakPoint from an existing SV 
         * @param sv pointer of SVRecord 
         * @param hdr bam header
         * @param cl predefined cc len
         */
        BreakPoint(const SVRecord* sv, const bam_hdr_t* hdr, int32_t cl){
            mSVStartBeg = sv->mSVStart;
            mSVStartEnd = sv->mSVStart;
            mSVEndBeg = sv->mSVEnd;
            mSVEndEnd = sv->mSVEnd;
            mSVStart = sv->mSVStart;
            mSVEnd = sv->mSVEnd;
            mPESupport = sv->mPESupport;
            mSRSupport = sv->mSRSupport;
            mSVT = sv->mSVT;
            mChr1 = sv->mChr1;
            mChr2 = sv->mChr2;
            mChr1Len = hdr->target_len[mChr1];
            mChr2Len = hdr->target_len[mChr2];
            mBoundary = cl;
            if(hdr) init(mChr1Len, mChr2Len);
        }
        /** BreakPoint default destructor */
        ~BreakPoint(){}

        /** operator to output an BreakPoint to ostream
         * @param os reference of ostream object
         * @param bp reference of BreakPoint object
         * @return reference of ostream object
         */
        inline friend std::ostream& operator<<(std::ostream& os, const BreakPoint& bp){
            os << "====================================\n";
            os << "SV Type: " << bp.mSVT << "\n";
            os << "SV larger chr: " << bp.mChr1 << "\n";
            os << "SV larger chr len: " << bp.mChr1Len << "\n";
            os << "SV little chr: " << bp.mChr2 << "\n";
            os << "SV little chr len: " << bp.mChr2Len << "\n";
            os << "SV leftmost starting position: " << bp.mSVStartBeg << "\n";
            os << "SV starting position: " << bp.mSVStart << "\n";
            os << "SV rightmost starting position: " << bp.mSVStartEnd << "\n";
            os << "SV leftmost ending position: " << bp.mSVEndBeg << "\n";
            os << "SV ending position: " << bp.mSVEnd << "\n";
            os << "SV rightmost ending position: " << bp.mSVEndEnd << "\n";
            os << "SV SR support number: " << bp.mSRSupport << "\n";
            os << "SV DP support number: " << bp.mPESupport << "\n";
            os << "SV SR consenus length: " << bp.mBoundary << "\n";
            os << "====================================\n";
            return os;
        }

    public:

        /** initialize an breakpoint to calculate its start/end boundary
         * @param largeChrLen mChr1 total reference length
         * @param smallChrLen mChr2 total reference length
         */
        void init(int32_t largeChrLen, int32_t smallChrLen);
        
        /** get sample reference sequence spanning an SV breakpoint[mSVStartBeg, mSVEndEnd]
         * @param smallChrSeq mChr2  total sequence
         * @param largeChrSeq mChr1 total sequence
         * @return contructed sample reference sequence spanning an SV breakpoint
         */
        std::string getSVRef(const char* smallChrSeq = NULL, const char* largeChrSeq = NULL);

        /** get larger chr seq of translocation SV
         * @param largeChrSeq mChr1 total sequence
         * @return large chr contributed reference sequence spanning an SV breakpoint
         */
        std::string getLargerRef(const char* largeChrSeq);

        /** get little chr seq of translocation SV
         * @param littleChrSeq mChr1 total sequence
         * @return little chr contributed reference sequence spanning an SV breakpoint
         */
        std::string getLittleRef(const char* littleChrSeq);
};

#endif
