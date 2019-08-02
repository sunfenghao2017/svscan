#ifndef ALIGNCFG_H
#define ALIGNCFG_H

#include <vector>
#include <memory>
#include <string>

/** class to store alignment score strategy */
class AlignConfig{
    public:
        int mMatch = 5;                     ///< score of same base alignment
        int mMisMatch = -4;                 ///< score of different base alignment
        int mGapOpen = -10;                 ///< score of gap open
        int mGapExt = -1;                   ///< score of gap extension
        int mInf = 1000000;                 ///< dummy variable for infinte
        bool mVerticalEndGapFree = false;   ///< vertical end gap costs 0 if true
        bool mHorizontalEndGapFree = false; ///< horizontal end gap costs 0 if true
    public:
        /** default constructor */
        AlignConfig(){}
        
        /** default destructor */
        ~AlignConfig(){}

        /** construct customized DnaScore object
         * @param match score of same base alignment
         * @param mismatch score of different base alignment
         * @param gapopen score of gap open
         * @param gapext score of gap open
         * @param vEndGapfree vertical end gap costs 0 if true
         * @param hEndGapFree horizontal end gap costs 0 if true
         */
        AlignConfig(int match, int mismatch, int gapopen, int gapext, bool vEndGapfree = false, bool hEndGapFree = false){
            mMatch = match;
            mMisMatch = mismatch;
            mGapOpen = gapopen;
            mGapExt = gapext;
            mInf = 1000000;
            mVerticalEndGapFree = vEndGapfree;
            mHorizontalEndGapFree = hEndGapFree;
        }

    public:
        /** get horizontal gap penalty accumulated
         * @param pos position at which to open/extend gap
         * @param end position at which the sequence ends
         * @param len gap length until pos
         */
        inline int horizontalGapSum(int pos, int end, int len){
            if(mHorizontalEndGapFree && (pos == 0 || pos == end)){
                return 0;
            }
            return mGapOpen + mGapExt * len;
        }

        /** get vertical gap penalty accumulated
         * @param pos position at which to open/extend gap
         * @param end position at which the sequence ends
         * @param len gap length until pos
         */
        inline int verticalGapSum(int pos, int end, int len){
            if(mVerticalEndGapFree && (pos == 0 || pos == end)){
                return 0;
            }
            return mGapOpen + mGapExt * len;
        }
        /** get horizontal gap penalty extend cost
         * @param pos position at which to extend gap
         * @param end position at which the sequence ends
         */
        inline int horizontalGapExtend(int pos, int end){
            if(mHorizontalEndGapFree && (pos == 0 || pos == end)){
                return 0;
            }
            return mGapExt;
        }

        /** get vertical gap penalty extend cost
         * @param pos position at which to extend gap
         * @param end position at which the sequence ends
         */
        inline int verticalGapExtend(int pos, int end){
            if(mVerticalEndGapFree && (pos == 0 || pos == end)){
                return 0;
            }
            return mGapExt;
        }
};


#endif
