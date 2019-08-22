#ifndef STATS_H
#define STATS_H

#include "util.h"
#include "svutil.h"
#include "options.h"
#include "svrecord.h"
#include "aligner.h"
#include "bamutil.h"
#include "matrix2d.h"
#include "aligncfg.h"
#include "srbamrecord.h"
#include "dpbamrecord.h"
#include "alndescriptor.h"
#include <unordered_map>
#include <htslib/sam.h>
#include <htslib/faidx.h>

/** class to store spanning PE supporting SV starting/ending positions */
struct SpanPoint{
    int32_t mBpPos = 0; ///< spanning PE supporting SV starting/ending positionss
    int32_t mSVT = 0;   ///< SV type
    int32_t mID = 0;    ///< SV ID, the index at which SV is stored in an vector

    /** constructor */
    SpanPoint(){}

    /** constructor
     * @param bpPos spanning PE supporting SV starting/ending positions
     */
    SpanPoint(int32_t bpPos) : mBpPos(bpPos) {}

    /** constructor
     * @param bpPos spanning PE supporting SV starting/ending positions
     * @param svt SV type
     * @param id SV ID, the index at which SV is stored in an vector
     */
    SpanPoint(int32_t bpPos, int32_t svt, int32_t id) : mBpPos(bpPos), mSVT(svt), mID(id) {};

    /** destructor */
    ~SpanPoint(){}

    /** operator to compare two SpanPoint object
     * @param other reference to another SpanPoint object
     * @return true if this < other
     */
    inline bool operator<(const SpanPoint& other) const {
        return mBpPos < other.mBpPos;
    }

    /** operator to output SpanPoint to ostream
     * @param os reference of ostream object
     * @param sp reference of SpanPoint object
     * @return reference of ostream object
     */
    inline friend std::ostream& operator<<(std::ostream& os, const SpanPoint& sp){
        os << "==================================\n";
        os << "Breakpoint position: " << sp.mBpPos << "\n";
        os << "SVT: " << sp.mSVT << "\n";
        os << "ID: " << sp.mID << "\n";
        os << "==================================\n";
        return os;
    }
};

/** class to store breakpoint region info */
struct BpRegion{
    int32_t mRegStart = 0; ///< region starting position on reference
    int32_t mRegEnd = 0;   ///< region ending position on reference
    int32_t mBpPos = 0;    ///< breakpoing position on reference
    int32_t mSVT = 0;      ///< SV type
    int32_t mID = 0;       ///< SV ID, the index at which SV is stored in an vector
    bool mIsSVEnd = false; ///< if true this region spanning SV ending position

    /** constructor */
    BpRegion(){}

    /** constructor
     * @param bpPos breakpoint position on ref
     */
    BpRegion(int32_t bpPos) : mBpPos(bpPos) {}

    /** destructor */
    ~BpRegion(){}

    /** operator to compare two BpRegion object
     * @param other reference to another BpRegion object
     * @return true if this < other
     */
    inline bool operator<(const BpRegion& other) const {
        return mBpPos < other.mBpPos;
    }

    /** operator to output BpRegion object to ostream
     * @param os reference of ostream object
     * @param br reference of BpRegion object
     * @return reference of ostream object
     */
    inline friend std::ostream& operator<<(std::ostream& os, const BpRegion& br){
        os << "======================================\n";
        os << "Region starting position: " << br.mRegStart << "\n";
        os << "Region ending position: " << br.mRegEnd << "\n";
        os << "Breakpoint position: " << br.mBpPos << "\n";
        os << "SVT: " << br.mSVT << "\n";
        os << "ID: " << br.mID << "\n";
        os << std::boolalpha << "Span SV ending position ? " << br.mIsSVEnd << "\n";
        os << "======================================\n";
        return os;
    }
};

/** Read coverage count of an SV event */
struct ReadCount{
    int32_t mLeftRC = 0;  ///< left side of SV starting position read count
    int32_t mRC = 0;      ///< read count between SV starting and ending
    int32_t mRightRC = 0; ///< right side of SV ending position read count

    /** constructor */
    ReadCount(){}

    /** destructor */
    ~ReadCount(){}

    /** operator to output ReadCount object to ostream
     * @param os reference of ostream object
     * @param rc reference of ReadCount object
     * @return reference of ostream object
     */
    inline friend std::ostream& operator<<(std::ostream& os, const ReadCount& rc){
        os << "=============================================================\n";
        os << "Reads count on left side of SV starting position: " << rc.mLeftRC << "\n";
        os << "Reads count between SV starting and ending position: " << rc.mRC << "\n";
        os << "Reads count on right side of SV ending position: " << rc.mRightRC << "\n";
        os << "=============================================================\n";
        return os;
    }
};

/** coverage range record of, each SV event needs 3 coverage range record\n
 * left control range, middle range, right control range\n
 * left control range, is the nearest left region which does not overlap with continuous N regions\n
 * middle region is just sv starting and ending range\n
 * right control range, is the nearest right region which does not overlap with continuous N regions
 * */
struct CovRecord{
    int32_t mStart = 0; ///< coverage starting position
    int32_t mEnd = 0;   ///< coverage ending position
    int32_t mID = 0;    ///< ID of this coverage range

    /** constructor */
    CovRecord(){}

    /** destructor */
    ~CovRecord(){}

    /** operator to compare two CovRecords
     * @param other reference of another CovRecord
     * @return true if this < other
     */
    inline bool operator<(const CovRecord& other) const {
        return mStart < other.mStart || 
               (mStart == other.mStart && mEnd < other.mEnd) || 
               (mStart == other.mStart && mEnd == other.mEnd && mID < other.mID);
    }

    /** operator to output CovRecord object to ostream
     * @param os reference of ostream object
     * @param cr reference of CovRecord object
     * @return reference of ostream object
     */
    inline friend std::ostream& operator<<(std::ostream& os, const CovRecord& cr){
        os << "======================================================\n";
        os << "Coverage counting region starting position: " << cr.mStart << "\n";
        os << "Coverage counting region ending posision: " << cr.mEnd << "\n";
        os << "SV ID of counting region: " << cr.mID << "\n";
        os << "======================================================\n";
        return os;
    }
};

/** Paired-end read which spanning SV breakpoint Stat object */
struct SpanningCount{
    int32_t mRefh1 = 0;            ///< count HP tag with value 1 supporting reference reads
    int32_t mRefh2 = 0;            ///< count HP tag with value not 1 supporting reference reads
    int32_t mAlth1 = 0;            ///< count HP tag with value 1 supporting SV reads
    int32_t mAlth2 = 0;            ///< count HP tag with value not 1 supporting SV reads
    std::vector<uint8_t> mRefQual; ///< Reads which are more likely to be reference mapping quality
    std::vector<uint8_t> mAltQual; ///< Reads which are more likely to be SV supporting mapping quality
    
    /** constructor */
    SpanningCount(){}

    /** destructor */
    ~SpanningCount(){}

    /** operator to output SpanningCount object to ostream
     * @param os reference of ostream object
     * @param sc reference of SpanningCount object
     * @return reference of ostream object
     */
    inline friend std::ostream& operator<<(std::ostream& os, const SpanningCount& sc){
        os << "=========================================================\n";
        os << "Read pairs which support REF haplotype: " << sc.mRefQual.size() << "\n";
        os << "Read pairs which support ALT haplotype: " << sc.mAltQual.size() << "\n";
        os << "HP == 1 read pair which support REF haplotype: " << sc.mRefh1 << "\n";
        os << "HP != 1 read pair which support REF haplotype: " << sc.mRefh2 << "\n";
        os << "HP == 1 read pair which support ALT haplotype: " << sc.mAlth1 << "\n";
        os << "HP != 1 read pair which support ALT haplotype: " << sc.mAlth2 << "\n";
        os << "=========================================================\n";
        return os;
    }
};

/** Single read which spanning SV breakpoint Stat object*/
struct JunctionCount{
    int32_t mRefh1 = 0;            ///< count HP tag with value 1 supporting reference reads
    int32_t mRefh2 = 0;            ///< count HP tag with value not 1 supporting reference reads
    int32_t mAlth1 = 0;            ///< count HP tag with value 1 supporting SV reads
    int32_t mAlth2 = 0;            ///< count HP tag with value not 1 supporting SV reads
    int32_t mFPIns = 0;            ///< false positive insertion supporting number
    std::vector<uint8_t> mRefQual; ///< Reads which are more likely to be reference mapping quality
    std::vector<uint8_t> mAltQual; ///< Reads which are more likely to be SV supporting mapping quality
    
    /** constructor */
    JunctionCount(){}

    /** destructor */
    ~JunctionCount(){} 
    
    /** operator to output JunctionCount object to ostream
     * @param os reference of ostream object
     * @param jc reference of JunctionCount object
     * @return reference of ostream object
     */
    inline friend std::ostream& operator<<(std::ostream& os, const JunctionCount& jc){
        os << "=========================================================\n";
        os << "Reas which support REF haplotype: " << jc.mRefQual.size() << "\n";
        os << "Reas which support ALT haplotype: " << jc.mAltQual.size() << "\n";
        os << "HP == 1 read which support REF haplotype: " << jc.mRefh1 << "\n";
        os << "HP != 1 read which support REF haplotype: " << jc.mRefh2 << "\n";
        os << "HP == 1 read which support ALT haplotype: " << jc.mAlth1 << "\n";
        os << "HP != 1 read which support ALT haplotype: " << jc.mAlth2 << "\n";
        os << "=========================================================\n";
        return os;
    }
};

/** class to store gene information of an SV */
class GeneInfo{
    public:
        std::string mGene1 = "-";         ///< gene name of breakpoint position on larger chrosome
        std::string mGene2 = "-";         ///< gene name of breakpoint position on small chrosome
        std::string mFuseGene = "-";      ///< fusion gene information(5'partner->3'partner)
        std::vector<std::string> mTrans1; ///< transcript names of breakpoint position on larger chrosome
        std::vector<std::string> mTrans2; ///< transcript names of breakpoint position on small chrosome

    public:
        /** default constructor */
        GeneInfo(){}

        /** default destructor */
        ~GeneInfo(){}
};

typedef std::vector<GeneInfo> GeneInfoList;
typedef std::vector<std::vector<BpRegion>> ContigBpRegions;
typedef std::vector<std::vector<SpanPoint>> ContigSpanPoints;

/** class to do coverage statistics of REF and ALT on contigs */
class Stats{
    public:
        Options* mOpt;                                     ///< pointer to Options
        int32_t mRefIdx;                                   ///< reference index
        std::vector<ReadCount> mReadCnts;                  ///< read count stat of each SV
        std::vector<JunctionCount> mJctCnts;               ///< Single read spanning SV breakpoint stats
        std::vector<SpanningCount> mSpnCnts;               ///< Paired-end read spanning SV breakpoint stats
        std::vector<std::pair<int32_t, int32_t>> mCovCnts; ///< base and fragment coverage count of each SV event
        std::vector<int32_t> mRefAlignedReadCount;         ///< REF like read count of each SV
        std::vector<int32_t> mRefAlignedSpanCount;         ///< REF like read pair count of each SV

    public:
        /** Stats constructor */
        Stats(){}

        /** Stats constructor
         * @param opt pointer to Options object
         * @param n total SVs to process
         * @param refidx reference index to compute statistics
         */
        Stats(Options* opt, int32_t n, int32_t refidx);

        /** Stats constructor
         * @param n total SV number
         */
        Stats(int32_t n);

        /** Stats destructor */
        ~Stats(){}

        /** create spaces
         * @param n SV total number
         */
        void init(int n);

        /** operator to output Stats result to ostream
         * @param os reference of ostream object
         * @param st reference of Stats object
         * @return reference of ostream object
         */
        inline friend std::ostream& operator<<(std::ostream& os, const Stats& st){
            for(uint32_t id = 0; id < st.mReadCnts.size(); ++id){
                os << "SV ID: " << id << "\n";
                os << "Read Count:\n" << st.mReadCnts[id];
                os << "Junction Read Count:\n" << st.mJctCnts[id];
                os << "Discordant Read Pair Count:\n" << st.mSpnCnts[id];
            }
            return os;
        }

        /** operator to output Stats result to ostream
         * @param os reference of ostream object
         * @param pst pointer to Stats object
         * @return reference of ostream object
         */
        inline friend std::ostream& operator<<(std::ostream& os, const Stats* pst){
            os << (*pst);
            return os;
        }

        /** gather coverage information of one contig
         * @param svs reference of SVSet(all SVs)
         * @param covRecs coverage records of 3-part of each SV events on each contig
         * @param bpRegs SV breakpoint regions on each contig
         * @param spPts SV DP read mapping position on each contig
         */
        void stat(const SVSet& svs, const std::vector<std::vector<CovRecord>>& covRecs,  const ContigBpRegions& bpRegs, const ContigSpanPoints& spPts);

        /** merge coverage information of all contigs
         * @param sts reference of list of Stats
         * @param n number of svs in total
         * @return merged stat info
         */
        static Stats* merge(const std::vector<Stats*>& sts, int32_t n);

        /** report BCF format report of all SVs
         * @param svs reference of SVSet
         */
        void reportBCF(const SVSet& svs);

        /** report TSV format report of all SVs
         * @param svs reference of SVSet
         * @param gl reference of GeneInfoList
         */
        void reportTSV(const SVSet& svs, const GeneInfoList& gl);

        /** get alignment quality of sequence against an read
         * @param alnResult align result(sequence vertical, read horizontal)
         * @param qual read quality
         * @return average mapping quality of matched bases
         */
        uint32_t getAlignmentQual(Matrix2D<char>* alnResult, const uint8_t* qual);

        /** test whether an bam record is met for the first time
         * @param b pointer to bam1_t struct
         * @param lastAlignedReads set of reads mapped at last position
         * @return true if b is met for the first time
         */
        inline static bool firstInPair(bam1_t* b, std::set<size_t>& lastAlignedReads){
            if(b->core.tid == b->core.mtid){
                return (b->core.pos < b->core.mpos) ||
                       (b->core.pos == b->core.mpos && lastAlignedReads.find(svutil::hashString(bam_get_qname(b))) == lastAlignedReads.end());
            }else return b->core.tid < b->core.mtid;
        }
};

#endif
