#ifndef SVRECORD_H
#define SVRECORD_H

#include <string>
#include <vector>
#include <cstdint>
#include "options.h"
#include "aligner.h"
#include "aligncfg.h"
#include "matrix2d.h"
#include "alndescriptor.h"
#include <htslib/faidx.h>
#include <htslib/sam.h>

/** class to store structural variant record */
class SVRecord{
    public:
        int32_t mChr1 = -1;          ///< chr on 5' end of SV | larger chr
        int32_t mSVStart = -1;       ///< 5' starting position of SV | starting position of SV on larger chr
        int32_t mChr2 = -1;          ///< chr on 3' end of SV | little chr
        int32_t mSVEnd = -1;         ///< 3' ending position of SV | ending position on little chr
        int32_t mCiPosLow = 0;       ///< smallest sv starting position is mSVStart + mCiPosLow(mCiPosLow is an negative value)
        int32_t mCiPosHigh = 0;      ///< largest sv starting position is mSVStart + mCiPosHigh(mCiPosHigh is an positive value)
        int32_t mCiEndLow = 0;       ///< smallest sv ending position is mSVEnd + mCiEndLow(mCiEndLow is an negative value)
        int32_t mCiEndHigh = 0;      ///< largest sv ending position is mSVEnd + mCiEndHigh(mCiEndHigh is an positive value)
        int32_t mPESupport = 0;      ///< number of Paired-end discordant reads supporting this SV
        int32_t mSRSupport = 0;      ///< number of Split-reads supporting this SV(one split read consists two part)
        int32_t mInsLen = 0;         ///< insertion size of this SV event due to Split-reads split alignment(absolute difference of seqpos of two part of SR)
        int32_t mHomLen = 0;         ///< total homology length of left/right part of consensus seq out of gap range with their gap elonged partner
        int32_t mSVT = -1;           ///< SV type[0-9]
        int32_t mID = -1;            ///< SV ID, is just the index at which this SVRecord is stored in the vector
        int32_t mSize = -1;          ///< SV size
        float mSRAlignQuality = 0;   ///< identity percent outside of largest inner gap in consensus SR sequence aligned with SV ref
        uint8_t mSRMapQuality = 0;   ///< median mapping quality of SR bam records that support this SV
        uint8_t mPEMapQuality = 0;   ///< median mapping quality of DP bam records that support this SV
        bool mPrecise = false;       ///< consensus sequence aligned with SV ref successfully and got a refined breakpoint evaluation
        std::string mAlleles = "";   ///< standard vcf format allele representation of SV
        std::string mConsensus = ""; ///< consensus sequence of SRs supporting this SV coming freom SRs MSA result
        std::string mSVRef = "";     ///< reference sequence of this SV constructed which spanning starting and ending positions
        int32_t mGapCoord[4] = {0};  ///< gap coordinates of split alignment of consensus sequence against reference[conGapBeg, conGapEnd, refGapBeg, refGapEnd]
        std::string mNameChr1 = "";  ///< name of chr on 5' end of SV | larger chr
        std::string mNameChr2 = "";  ///< name of chr on 3' end of SV | little chr
        std::string mProbeBegC = ""; ///< an consensus sequence segment spanning the SV starting position
        std::string mProbeBegR = ""; ///< an reference sequence segment spanning the SV starting position
        std::string mProbeEndC = ""; ///< an consensus sequence segment spanning the SV ending position
        std::string mProbeEndR = ""; ///< an reference sequence segment spanning the SV ending position

    public:
        /** SVRecord constructor */
        SVRecord() {}

        /** SVRecord destructor */
        ~SVRecord() {}

    public:
        /** operator to output SVRecord to ostream
         * @param os reference of ostream object
         * @param sv reference of SVRecord object
         * @return reference of ostream object
         */
        inline friend std::ostream& operator<<(std::ostream& os, const SVRecord& sv){
            os << "======================================================================================================\n";
            os << "SV type: " << sv.mSVT << "\n";
            os << "SV ID: " << sv.mID << "\n";
            os << "Reference ID on 5' end of SV: " << sv.mChr1 << "(" << sv.mNameChr1 << ")\n";
            os << "Starting position on reference of 5' end of SV: " << sv.mSVStart << "\n";
            os << "Reference ID on 3' end of SV: " << sv.mChr2 << "(" << sv.mNameChr2 << ")\n";
            os << "Ending position on reference of 3' end of SV: " << sv.mSVEnd << "\n";
            os << "Negative offset of starting position on reference of 5' end of SV: " << sv.mCiPosLow << "\n";
            os << "Positive offset of starting position on reference of 3' end of SV: " << sv.mCiPosHigh << "\n";
            os << "Negative offset of ending position on reference of 5' end of SV: " << sv.mCiEndLow << "\n";
            os << "Positive offset of ending position on reference of 3' end of SV: " << sv.mCiEndHigh << "\n";
            os << "Number of discordant paired-end reads supporting this SV: " << sv.mPESupport << "\n";
            os << "Number of split reads supporting this SV: " << sv.mSRSupport << "\n";
            os << "Insertion size of this SV contributed by split read: " << sv.mInsLen << "\n";
            os << "Total homology length between both end of consensus read vs the reference sequence: " << sv.mHomLen << "\n";
            os << "Identity percentage of consensus split read against reference outside of inner longest gap: " << sv.mSRAlignQuality << "\n";
            os << "Median mapping quality of all split read alignment record which support this SV: " << sv.mSRMapQuality << "\n";
            os << "Median mapping quality of all discordant paired-end reads alignment record which support this SV: " << sv.mPEMapQuality << "\n";
            os << std::boolalpha << "Consensus split read split aligned against reference got a refined breakpoint: " << sv.mPrecise << "\n";
            os << "Allele of this SV event: " << sv.mAlleles << "\n";
            os << "Consensus sequence of split reads supporting this SV: " << sv.mConsensus << "\n";
            os << "Constructed reference sequence of this SV: " << sv.mSVRef << "\n";
            os << "Consensus sequence segment spanning the SV starting position: " << sv.mProbeBegC << "\n";
            os << "Consensus sequence segment spanning the SV ending position: " << sv.mProbeEndC << "\n";
            os << "Reference sequence segment spanning the SV starting position: " << sv.mProbeBegR << "\n";
            os << "Reference sequence segment spanning the SV ending position: " << sv.mProbeEndR << "\n";
            os << "======================================================================================================\n";
            return os;
        }

        /** operator to compare two SVRecords
         * @param other reference to another SVRecord
         * @return true if this < other
         */
        inline bool operator<(const SVRecord& other) const {
            return mChr1 < other.mChr1 || (mChr1 == other.mChr1 && mSVStart < other.mSVStart) ||
                   (mChr1 == other.mChr1 && mSVStart == other.mSVStart && mSVEnd < other.mSVEnd) ||
                   (mChr1 == other.mChr1 && mSVStart == other.mSVStart && mSVEnd == other.mSVEnd && mPESupport > other.mPESupport);
        }

        /** align consensus SR seq to constructed SV ref seq to refine breakpoint coordinate
         * @param opt pointer to Options object
         * @param hdr bam header
         * @param chr1Seq little chr sequence
         * @param chr2Seq larger chr sequence
         * @return true if breakpoint refined
         */
        bool refineSRBp(const Options* opt, const bam_hdr_t* hdr, const char* chr1Seq, const char* chr2Seq = NULL);
        
        /** align consensus SR seq to constructed SV ref seq by split alignment strategy
         * @param alnResult to storealignment result
         * @return true if a fine split alignment result found
         */
        bool consensusRefAlign(Matrix2D<char>* alnResult);
        
        /** find new gap range from split alignment result
         * @param alnResult consensus SR and constructed SV ref split alignment result
         * @param ad reference of AlignDescriptor to store refined align info
         * @return true if a better split found
         */
        bool findSplit(Matrix2D<char>* alnResult, AlignDescriptor& ad);
        
        /** find homology length of part of SR consensus out of inner largest gap against elonged SV ref in gap
         * @param ad reference of AlignDescriptor to store result
         */
        void findHomology(AlignDescriptor& ad);
        
        /** transform refined gap coordinates to genome coordinate
         * @param bp BreakPoint object storing initialized breakpoing coordinates(which used to extrace ref seq of SV)
         * @param ad reference of AlignDescriptor storing refined gap alignment info
         * @param finalGapStart to store final gap starting position
         * @param finalGapEnd to store final gap ending position
         */
        template<typename TBreakPoint>
        bool coordTransform(TBreakPoint& bp, AlignDescriptor& ad, int32_t& finalGapStart, int32_t& finalGapEnd);

        /** check if new continous gaps opened in middle of SR consensus and sv ref is valid
         * @param refGap new sv ref gap length
         * @param oldRefGap old sv ref gap length
         * @param varGap new SR consensus gap length
         * @param oldVarGap old SR consensus gap length
         */
        inline bool checkSVGap(int32_t refGap, int32_t oldRefGap, int32_t varGap, int32_t oldVarGap){
            if(mSVT == 4) return varGap > oldVarGap;
            return refGap > oldRefGap;
        }

        /** parse SVT marker to string marker
         * @return str representation of SVT
         */
        inline std::string addID(){
            if(mSVT == 0 || mSVT == 1) return "INV";
            else if(mSVT == 2) return "DEL";
            else if(mSVT == 3) return "DUP";
            else if(mSVT == 4) return "INS";
            else return "BND";
        };

        /** add allele information of this SV
         * @param ref reference sequence of mChr1 in [mSVStart-1, mSVStart]
         */
        inline void addAlleles(){
            std::string ref = mSVRef;
            if(mPrecise) ref = mSVRef.substr(mGapCoord[2] - 1, 1);
            if(mSVT >= 5){
                if(mSVT == 5){// 5to5 translocation
                    mAlleles = ref + "," + ref + "]" + mNameChr2 + ":" + std::to_string(mSVEnd) + "]";
                }else if(mSVT == 6){// 3to3 translocation
                    mAlleles = ref + "," + "[" + mNameChr2 + ":" + std::to_string(mSVEnd) + "[" + ref;
                }else if(mSVT == 7){// 5to3 translocation
                    mAlleles = ref + "," + ref + "[" + mNameChr2 + ":" + std::to_string(mSVEnd) + "[";
                }else if(mSVT == 8){// 3to5 translocation
                    mAlleles = ref + "," + "]" + mNameChr2 + ":" + std::to_string(mSVEnd) + "]" + ref;
                }else{
                    mAlleles = ref + ",<" + addID() + ">";
                }
            }else mAlleles = ref + ",<" + addID() + ">";
        }
};

/** type to store a list of structural variant record */
typedef std::vector<SVRecord> SVSet; ///< list of SV

inline std::ostream& operator<<(std::ostream& os, const SVSet& svs){
    for(uint32_t id = 0; id < svs.size(); ++id) os<< svs[id];
    return os;
}

/** merge SVSet supported by SR and PE, refine PE SV with SR SV and keep better SR SV if dup exists
 * @param sr SVSet supported by SR
 * @param pe SVSet supported by PE
 * @param peMergeSearchWindow the SV starting pos range to search for each SV in sr for PE support
 * @param srDupSearchWindow the SV starting pos range to search for an non-PE support SR SV for dup
 */
void mergeAndSortSVSet(SVSet& sr, SVSet& pe, int32_t peMergeSearchWindow = 500, int32_t srDupSearchWindow = 10);

/** get reference of  SV supported by PE
 * @param pe SVSet supported by PE
 * @param opt pointer to Options object
 * @param boundary max length offset at each breakpoint to fetch reference
 */
void getDPSVRef(SVSet& pe, Options* opt);

#endif
