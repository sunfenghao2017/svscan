#include "breakpoint.h"

void BreakPoint::init(int32_t largeChrLen, int32_t smallChrLen){
    if(mSVT >= 4){// insetion and translocation
        mSVStartBeg = std::max(0, mSVStart - mBoundary);
        mSVStartEnd = std::min(largeChrLen, mSVStart + mBoundary);
        mSVEndBeg = std::max(0, mSVEnd - mBoundary);
        mSVEndEnd = std::min(smallChrLen, mSVEnd + mBoundary);
    }else{
        mSVStartBeg = std::max(0, mSVStart - mBoundary);
        mSVStartEnd = std::min(mSVStart + mBoundary, (mSVStart + mSVEnd)/2);
        mSVEndBeg = std::max(mSVEnd - mBoundary, (mSVStart + mSVEnd)/2 + 1);
        mSVEndEnd = std::min(smallChrLen, mSVEnd + mBoundary);
    }
}

std::string BreakPoint::getSVRef(const char* smallChrSeq, const char* largeChrSeq){
    std::string chr2Part;
    // Get chr2 seq for tranclocation
    if(mSVT == 5){
        chr2Part = std::string(smallChrSeq + mSVEndBeg, smallChrSeq + mSVEndEnd);
        util::reverseComplement(chr2Part);
    }else if(mSVT > 5){
        chr2Part = std::string(smallChrSeq + mSVEndBeg, smallChrSeq + mSVEndEnd);
        util::str2upper(chr2Part);
    }
    // Get full ref seq for translocation
    if(mSVT >= 5){
        if(mSVT == 5 || mSVT == 7){// 5to5 || 5to3
            std::string chr1Part = std::string(largeChrSeq + mSVStartBeg, largeChrSeq + mSVStartEnd);
            util::str2upper(chr1Part);
            return chr1Part + chr2Part;
        }else if(mSVT == 6){// 3to3
            std::string chr1Part = std::string(largeChrSeq + mSVStartBeg, largeChrSeq + mSVStartEnd);
            util::reverseComplement(chr1Part);
            return chr1Part + chr2Part;
        }else{// 3to5
            std::string chr1Part = std::string(largeChrSeq + mSVStartBeg, largeChrSeq + mSVStartEnd);
            util::str2upper(chr1Part);
            return chr2Part + chr1Part;
        }
    }
    // Get full ref seq for SV on same chr
    if(mSVT == 0){// 5to5 left breakpoint of inversion
        std::string refEnd = std::string(smallChrSeq + mSVEndBeg, smallChrSeq + mSVEndEnd);
        util::reverseComplement(refEnd);
        std::string refBeg = std::string(smallChrSeq + mSVStartBeg, smallChrSeq + mSVStartEnd);
        util::str2upper(refBeg);
        return refBeg + refEnd;
    }
    if(mSVT == 1){// 3to3 right breakpoing of inversion
        std::string refBeg = std::string(smallChrSeq + mSVStartBeg, smallChrSeq + mSVStartEnd);
        util::reverseComplement(refBeg);
        std::string refEnd = std::string(smallChrSeq + mSVEndBeg, smallChrSeq + mSVEndEnd);
        util::str2upper(refEnd);
        return refBeg + refEnd;
    }
    if(mSVT == 2){// 5to3 Deletion
        std::string refBeg = std::string(smallChrSeq + mSVStartBeg, smallChrSeq + mSVStartEnd);
        std::string refEnd = std::string(smallChrSeq + mSVEndBeg, smallChrSeq + mSVEndEnd);
        util::str2upper(refBeg);
        util::str2upper(refEnd);
        return refBeg + refEnd;
    }
    if(mSVT == 3){// 3to5 Duplication
        std::string refBeg = std::string(smallChrSeq + mSVStartBeg, smallChrSeq + mSVStartEnd);
        std::string refEnd = std::string(smallChrSeq + mSVEndBeg, smallChrSeq + mSVEndEnd);
        util::str2upper(refBeg);
        util::str2upper(refEnd);
        return refEnd + refBeg;
    }
    if(mSVT == 4){// Insertion
        std::string refSeq = std::string(smallChrSeq + mSVStartBeg, smallChrSeq + mSVEndEnd);
        util::str2upper(refSeq);
        return refSeq;
    }
    return "";
}
