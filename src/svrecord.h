#ifndef SVRECORD_H
#define SVRECORD_H

#include <string>
#include <vector>
#include <cstdint>
#include "bamutil.h"
#include "options.h"
#include "aligner.h"
#include "aligncfg.h"
#include "matrix2d.h"
#include "realnfilter.h"
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
        int32_t mAlnInsLen = 0;      ///< insertion size of this SV event due to Split-reads split alignment(absolute difference of seqpos of two part of SR)
        std::string mBpInsSeq = "";  ///< insertion sequence after breakpoint position
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
        std::string mInsSeq = "";    ///< insertion sequence of insertion event
        std::string mTraChr1Seq = "";///< larger chr reference sequence of translocation 
        std::string mTraChr2Seq = "";///< little chr reference sequence of translocation
        bool mMerged = false;        ///< this SV has been merged if true
        bool mFromOneSR = false;     ///< this SV comes from one seed SR
        bool mPassRealn = true;      ///< this SV has passed realignment by online bwa

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
            os << "Beg position on reference of 5' end of SV: " << sv.mSVStart << "\n";
            os << "Reference ID on 3' end of SV: " << sv.mChr2 << "(" << sv.mNameChr2 << ")\n";
            os << "Ending position on reference of 3' end of SV: " << sv.mSVEnd << "\n";
            os << "Negative offset of starting position on reference of 5' end of SV: " << sv.mCiPosLow << "\n";
            os << "Positive offset of starting position on reference of 3' end of SV: " << sv.mCiPosHigh << "\n";
            os << "Negative offset of ending position on reference of 5' end of SV: " << sv.mCiEndLow << "\n";
            os << "Positive offset of ending position on reference of 3' end of SV: " << sv.mCiEndHigh << "\n";
            os << "Number of discordant paired-end reads supporting this SV: " << sv.mPESupport << "\n";
            os << "Number of split reads supporting this SV: " << sv.mSRSupport << "\n";
            os << "Insertion size of this SV contributed by split read: " << sv.mAlnInsLen << "\n";
            os << "Total homology length between both end of consensus read vs the reference sequence: " << sv.mHomLen << "\n";
            os << "Identity percentage of consensus split read against reference outside of inner longest gap: " << sv.mSRAlignQuality << "\n";
            os << "Median mapping quality of all split read alignment record which support this SV: " << sv.mSRMapQuality << "\n";
            os << "Median mapping quality of all discordant paired-end reads alignment record which support this SV: " << sv.mPEMapQuality << "\n";
            os << "Consensus split read split aligned against reference got a refined breakpoint: " << std::boolalpha << sv.mPrecise << "\n";
            os << "Merged by other SV event: " << std::boolalpha << sv.mMerged << "\n";
            os << "Allele of this SV event: " << sv.mAlleles << "\n";
            os << "Consensus sequence of split reads supporting this SV: " << sv.mConsensus << "\n";
            os << "Constructed reference sequence of this SV: " << sv.mSVRef << "\n";
            os << "Consensus sequence segment spanning the SV starting position: " << sv.mProbeBegC << "\n";
            os << "Consensus sequence segment spanning the SV ending position: " << sv.mProbeEndC << "\n";
            os << "Reference sequence segment spanning the SV starting position: " << sv.mProbeBegR << "\n";
            os << "Reference sequence segment spanning the SV ending position: " << sv.mProbeEndR << "\n";
            os << "Translocation chr1Seq: " << sv.mTraChr1Seq << "\n";
            os << "Translocation chr2Seq: " << sv.mTraChr2Seq << "\n";
            if(sv.mSVT == 4) os << "Inserted sequence: " << sv.mInsSeq << "\n";
            if(sv.mBpInsSeq.length() > 0) os << "Sequence inserted after break point: " << sv.mBpInsSeq << "\n";
            os << "======================================================================================================\n";
            return os;
        }

        /** convert an SVRecord object to string
         * @return str representation of SVRecord
         */
        inline std::string toStr() const {
            std::stringstream ss;
            ss << "\n======================================================================================================\n";
            ss << "SV type: " << mSVT << "\n";
            ss << "SV ID: " << mID << "\n";
            ss << "Reference ID on 5' end of SV: " << mChr1 << "(" << mNameChr1 << ")\n";
            ss << "Beg pssition on reference of 5' end of SV: " << mSVStart << "\n";
            ss << "Reference ID on 3' end of SV: " << mChr2 << "(" << mNameChr2 << ")\n";
            ss << "Ending pssition on reference of 3' end of SV: " << mSVEnd << "\n";
            ss << "Negative offset of starting pssition on reference of 5' end of SV: " << mCiPosLow << "\n";
            ss << "Positive offset of starting pssition on reference of 3' end of SV: " << mCiPosHigh << "\n";
            ss << "Negative offset of ending pssition on reference of 5' end of SV: " << mCiEndLow << "\n";
            ss << "Pssitive offset of ending pssition on reference of 3' end of SV: " << mCiEndHigh << "\n";
            ss << "Number of discordant paired-end reads supporting this SV: " << mPESupport << "\n";
            ss << "Number of split reads supporting this SV: " << mSRSupport << "\n";
            ss << "Insertion size of this SV contributed by split read: " << mAlnInsLen << "\n";
            ss << "Total homology length between both end of consensus read vs the reference sequence: " << mHomLen << "\n";
            ss << "Identity percentage of consensus split read against reference outside of inner longest gap: " << mSRAlignQuality << "\n";
            ss << "Median mapping quality of all split read alignment record which support this SV: " << mSRMapQuality << "\n";
            ss << "Median mapping quality of all discordant paired-end reads alignment record which support this SV: " << mPEMapQuality << "\n";
            ss << "Consensus split read split aligned against reference got a refined breakpoint: " << std::boolalpha << mPrecise << "\n";
            ss << "Merged by other SV event: " << std::boolalpha << mMerged << "\n";
            ss << "Allele of this SV event: " << mAlleles << "\n";
            ss << "Consensus sequence of split reads supporting this SV: " << mConsensus << "\n";
            ss << "Constructed reference sequence of this SV: " << mSVRef << "\n";
            ss << "Consensus sequence segment spanning the SV starting pssition: " << mProbeBegC << "\n";
            ss << "Consensus sequence segment spanning the SV ending pssition: " << mProbeEndC << "\n";
            ss << "Reference sequence segment spanning the SV starting pssition: " << mProbeBegR << "\n";
            ss << "Reference sequence segment spanning the SV ending pssition: " << mProbeEndR << "\n";
            ss << "Translocation chr1Seq: " << mTraChr1Seq << "\n";
            ss << "Translocation chr2Seq: " << mTraChr2Seq << "\n";
            if(mSVT == 4) ss << "Inserted sequence: " << mInsSeq << "\n";
            if(mBpInsSeq.length() > 0) ss << "Sequence inserted after break point: " << mBpInsSeq << "\n";
            ss << "======================================================================================================\n";
            return ss.str();
        }


        /** operator to compare two SVRecords
         * @param other reference to another SVRecord
         * @return true if this < other
         */
        inline bool operator<(const SVRecord& other) const {
            return (mSVT < other.mSVT) ||
                   (mSVT == other.mSVT && mChr1 < other.mChr1) || 
                   (mSVT == other.mSVT && mChr1 == other.mChr1 && mChr2 < other.mChr2) ||
                   (mSVT == other.mSVT && mChr1 == other.mChr1 && mChr2 == other.mChr2 && mSVStart < other.mSVStart) ||
                   (mSVT == other.mSVT && mChr1 == other.mChr1 && mChr2 == other.mChr2 && mSVStart == other.mSVStart && mSVEnd < other.mSVEnd) || 
                   (mSVT == other.mSVT && mChr1 == other.mChr1 && mChr2 == other.mChr2 && mSVStart == other.mSVStart && mSVEnd == other.mSVEnd && mSRSupport < other.mSRSupport) ||
                   (mSVT == other.mSVT && mChr1 == other.mChr1 && mChr2 == other.mChr2 && mSVStart == other.mSVStart && mSVEnd == other.mSVEnd && mSRSupport == other.mSRSupport && mPESupport < other.mPESupport);
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

         /** get read sequence with inserted sequence extarcted away
          * @param b pointer to bam1_t struct
          * @param rseq read seq without inserted sequences
          * @param iseq inserted sequence
          * @param inslen inserted length after clip pos
          * @param maxReadSep max inserted length allowed
          */
        inline void getSCIns(const bam1_t* b, std::string& rseq, std::string& iseq, int32_t inslen, int32_t maxReadSep){
            if(mSVT == 4 || inslen <= maxReadSep ){
                rseq = bamutil::getSeq(b);
                return;
            }
            std::string wseq = bamutil::getSeq(b);
            std::pair<int, int> sclen = bamutil::getSoftClipLength(b);
            if(sclen.first){
                if(sclen.first <= inslen) return;
                rseq = wseq.substr(0, sclen.first - inslen);
                rseq.append(wseq.substr(sclen.first));
                iseq = wseq.substr(sclen.first - inslen, inslen);
            }else if(sclen.second){
                if(sclen.second <= inslen) return;
                rseq = wseq.substr(0, b->core.l_qseq - sclen.second);
                rseq.append(wseq.substr(b->core.l_qseq - sclen.second + inslen));
                iseq = wseq.substr(b->core.l_qseq - sclen.second, inslen);
            }
        }

        /** add allele information of this SV */
        inline void addAlleles(){
            if(mSVT >= 5){
                if(mSVT == 5){// 5to5 translocation
                    mAlleles = mSVRef + "," + mSVRef + "]" + mNameChr2 + ":" + std::to_string(mSVEnd) + "]";
                }else if(mSVT == 6){// 3to3 translocation
                    mAlleles = mSVRef + "," + "[" + mNameChr2 + ":" + std::to_string(mSVEnd) + "[" + mSVRef;
                }else if(mSVT == 7){// 5to3 translocation
                    mAlleles = mSVRef + "," + mSVRef + "[" + mNameChr2 + ":" + std::to_string(mSVEnd) + "[";
                }else if(mSVT == 8){// 3to5 translocation
                    mAlleles = mSVRef + "," + "]" + mNameChr2 + ":" + std::to_string(mSVEnd) + "]" + mSVRef;
                }else{
                    mAlleles = mSVRef + ",<" + addID() + ">";
                }
            }else mAlleles = mSVRef + ",<" + addID() + ">";
        }
        
        /** merge two part of ref of translocation */
        inline void mergeRef(){
            if(mSVT == 8){
                mSVRef = mTraChr2Seq + mTraChr1Seq;
            }else{
                mSVRef = mTraChr1Seq + mTraChr2Seq;
            }
        }

        /** add reference probes of translocation SV */
        inline void addRefProbe(const Options* opt){
            if(mSVT == 8){
                int32_t beg = std::max(0, mGapCoord[3] - opt->filterOpt->mMinFlankSize);
                int32_t len = std::min((int32_t)mSVRef.length() - beg, 2 * opt->filterOpt->mMinFlankSize);
                mProbeBegR = mSVRef.substr(beg, len);
                beg = std::max(0, mGapCoord[2] - opt->filterOpt->mMinFlankSize);
                len = std::min((int32_t)mTraChr2Seq.length(), 2 * opt->filterOpt->mMinFlankSize);
                mProbeEndR = mSVRef.substr(beg, len);
            }else{
                int32_t beg = std::max(0, mGapCoord[2] - opt->filterOpt->mMinFlankSize);
                int32_t len = std::min((int32_t)mTraChr1Seq.length() - beg, 2 * opt->filterOpt->mMinFlankSize);
                mProbeBegR = mSVRef.substr(beg, len);
                beg = std::max(0, mGapCoord[3] - opt->filterOpt->mMinFlankSize);
                len = std::min((int32_t)mSVRef.length(), 2 * opt->filterOpt->mMinFlankSize);
                mProbeEndR = mSVRef.substr(beg, len);
            }
            if(mSVT == 6){
                util::reverseComplement(mProbeBegR);
            }
            if(mSVT == 5){
                util::reverseComplement(mProbeEndR);
            }
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
 * @param svs SVSet to store merged SVs
 * @param opt pointer to Options
 */
void mergeAndSortSVSet(SVSet& sr, SVSet& pe, SVSet& svs, Options* opt);

/** merge SVSet supported by SR only
 * @param sr SVSet supported by SR
 * @param msr SVSet merged
 * @param opt pointer to Options
 */
void mergeSRSVs(SVSet& sr, SVSet& msr, Options* opt);

/** merge SVSet supported by DP only
 * @param dp SVSet supported by SR
 * @param mdp SVSet merged
 * @param opt pointer to Options
 */
void mergeDPSVs(SVSet& dp, SVSet& mdp, Options* opt);

/** get reference of SV supported by PE
 * @param pe SVSet supported by PE
 * @param opt pointer to Options object
 * @param boundary max length offset at each breakpoint to fetch reference
 */
void getDPSVRef(SVSet& pe, Options* opt);

#endif
