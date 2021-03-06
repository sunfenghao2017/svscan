#ifndef STATS_H
#define STATS_H

#include <util.h>
#include "trsrec.h"
#include "svutil.h"
#include "fuserec.h"
#include "options.h"
#include "svrecord.h"
#include "aligner.h"
#include <bamutil.h>
#include "matrix2d.h"
#include "aligncfg.h"
#include "srbamrecord.h"
#include "dpbamrecord.h"
#include "alndescriptor.h"
#include <unordered_map>
#include <htslib/sam.h>
#include <htslib/faidx.h>

/** class to store reads count in each contig */
struct RegItemCnt{
    int64_t mCount = 0;       ///< read counts
    int32_t mTid = -1;        ///< contig index
    int32_t mBeg = -1;        ///< contig beg
    int32_t mEnd = -1;        ///< contig end
    bool mInterleved = false; ///< if true this region is contiguous to previous one

    /** operator to compare two contig read counts
     * @param other reference of RegItemCnt
     * @return true it this one has less read counts
     */
    inline bool operator<(const RegItemCnt& other) const{
        return mCount > other.mCount;
    }

    /** operator to output a RegItemCnt to ostream
     * @param os reference of ostream
     * @param ric reference of RegItemCnt 
     * @return reference of ostream
     */
    inline friend std::ostream& operator<<(std::ostream& os, const RegItemCnt& ric){
        os << ric.mTid << "\t[" << ric.mBeg << "," << ric.mEnd << "] " << ric.mCount << "," << std::boolalpha << ric.mInterleved << std::endl;
        return os;
    }
};

/** class to store spanning PE supporting SV starting/ending positions */
struct SpanPoint{
    int32_t mBpPos = 0;    ///< spanning PE supporting SV starting/ending positionss
    int32_t mSVT = 0;      ///< SV type
    int32_t mID = 0;       ///< SV ID, the index at which SV is stored in an vector
    bool mIsSVEnd = false; ///< end of sv pos if true

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
    SpanPoint(int32_t bpPos, int32_t svt, int32_t id, bool svEnd) : mBpPos(bpPos), mSVT(svt), mID(id), mIsSVEnd(svEnd) {};

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
    int32_t mRefCntBeg = 0;                 ///< number of reads which are more likely to be reference on sv starting pos
    int32_t mRefCntEnd = 0;                 ///< number of reads which are more likely to be reference on sv ending pos
    int32_t mAltCnt = 0;                    ///< number of reads which are more likely to be SV supporting on sv starting pos
    std::map<uint8_t, int32_t> mRefQualBeg; ///< mapQ of reads which are more likely to be reference on sv starting pos
    std::map<uint8_t, int32_t> mRefQualEnd; ///< mapQ of reads which are more likely to be reference on sv ending pos
    std::map<uint8_t, int32_t> mAltQual;    ///< mapQ of reads which are more likely to be SV supporting
    
    /** constructor */
    SpanningCount(){}

    /** destructor */
    ~SpanningCount(){}
    
    /** get reference depth of one sv event
     * @return ref depth of sv
     */
    inline int32_t getRefDep(){
        if(mRefCntBeg > mRefCntEnd) return mRefCntBeg;
        else return mRefCntEnd;
    }
    
    /** get alt depth of one sv event
     * @return alt depth of sv
     */
    inline int32_t getAltDep(){
        return mAltCnt;
    }

    /** operator to output SpanningCount object to ostream
     * @param os reference of ostream object
     * @param sc reference of SpanningCount object
     * @return reference of ostream object
     */
    inline friend std::ostream& operator<<(std::ostream& os, const SpanningCount& sc){
        os << "=========================================================\n";
        os << "Read pairs which support REF haplotype on SV starting pos: " << sc.mRefCntBeg << "\n";
        os << "Read pairs which support REF haplotype on SV ending pos: " << sc.mRefCntEnd << "\n";
        os << "Read pairs which support ALT haplotype: " << sc.mAltCnt << "\n";
        os << "=========================================================\n";
        return os;
    }
};

/** Single read which spanning SV breakpoint Stat object*/
struct JunctionCount{
    int32_t mFPIns = 0;                     ///< false positive insertion supporting number
    int32_t mRefCntBeg = 0;                 ///< number of reads which are more likely to be reference on sv starting pos
    int32_t mRefCntEnd = 0;                 ///< number of reads which are more likely to be reference on sv ending pos
    int32_t mAltCnt = 0;                    ///< number of reads which are more likely to be SV supporting on sv starting pos
    std::map<uint8_t, int32_t> mRefQualBeg; ///< mapQ of reads which are more likely to be reference on sv starting pos
    std::map<uint8_t, int32_t> mRefQualEnd; ///< mapQ of reads which are more likely to be reference on sv ending pos
    std::map<uint8_t, int32_t> mAltQual;    ///< mapQ of reads which are more likely to be SV supporting
    
    /** constructor */
    JunctionCount(){}

    /** destructor */
    ~JunctionCount(){} 

    /** get reference depth of one sv event
     * @return ref depth of sv
     */
    inline int32_t getRefDep(){
        if(mRefCntBeg > mRefCntEnd) return mRefCntBeg;
        else return mRefCntEnd;
    }
    
    /** get alt depth of one sv event
     * @return alt depth of sv
     */
    inline int32_t getAltDep(){
        return mAltCnt;
    }
    
    /** operator to output JunctionCount object to ostream
     * @param os reference of ostream object
     * @param jc reference of JunctionCount object
     * @return reference of ostream object
     */
    inline friend std::ostream& operator<<(std::ostream& os, const JunctionCount& jc){
        os << "=========================================================\n";
        os << "Split reads which support REF haplotype on SV starting pos: " << jc.mRefCntBeg << "\n";
        os << "Split reads which support REF haplotype on SV ending pos: " << jc.mRefCntEnd << "\n";
        os << "Split reads which support ALT haplotype: " << jc.mAltCnt << "\n";
        os << "=========================================================\n";
        return os;
    }
};

/** class to store a list of FuseGene */
typedef std::vector<FuseGene> FuseGeneList;

/** class to store gene information of an SV */
class GeneInfo{
    public:
        TrsRecList mGene1;         ///< gene transcript records of breakpoint on larger chrosome
        TrsRecList mGene2;         ///< gene transcript records of breakpoint on little chrosome
        FuseGeneList mFuseGene;    ///< fusion gene information
        int32_t mPos1 = 0;         ///< position on genome (RNA sv only)
        int32_t mPos2 = 0;         ///< position on genome (RNA sv only)
        std::string mChr1 = ".";   ///< chr of gene1
        std::string mChr2 = ".";   ///< chr of gene2
        std::string mFsCigar = ""; ///< fsCigar (RNA sv only)

    public:
        /** default constructor */
        GeneInfo(){}

        /** default destructor */
        ~GeneInfo(){}

        /** operator output GeneInfo to stream
         * @param os reference of ostream
         * @param gi reference of GeneInfo
         * @return reference of ostream
         */
        inline friend std::ostream& operator<<(std::ostream& os, GeneInfo& gi){
            os << "Gene1: ";
            for(uint32_t i = 0; i < gi.mGene1.size(); ++i){
                os << gi.mGene1[i].toStr();
                if(i != gi.mGene1.size() - 1) os << ";";
            }
            os << "\n";
            os << "Gene2: ";
            for(uint32_t i = 0; i < gi.mGene2.size(); ++i){
                os << gi.mGene2[i].toStr();
                if(i != gi.mGene2.size() - 1) os << ";";
            }
            os << "\n";
            os << "FuseGene: ";
            for(uint32_t i = 0; i < gi.mFuseGene.size(); ++i){
                os << gi.mFuseGene[i];
                if(i != gi.mFuseGene.size() - 1) os << ",";
            }
            os << "\n";
            os << "Pos1: " << gi.mPos1 << "\n";
            os << "Pos2: " << gi.mPos2 << "\n";
            os << "Chr1: " << gi.mChr1 << "\n";
            os << "Chr2: " << gi.mChr2 << "\n";
            return os;
        }

        /** get str rep of FuseGene
         * @return string rep of FuseGene
         */
        inline std::string getFuseGene(){
            std::vector<std::string> fstr;
            for(uint32_t i = 0; i < mFuseGene.size(); ++i){
                fstr.push_back(mFuseGene[i].toStr());
            }
            return util::join(fstr, ";");
        }

        /** get fsmask of FuseGene
         * @return fsmask of FuseGene
         */
        inline std::string getFsMask(){
            std::vector<std::string> fmsk;
            for(uint32_t i = 0; i < mFuseGene.size(); ++i){
                fmsk.push_back(std::to_string(mFuseGene[i].status));
            }
            return util::join(fmsk, ";");
        }

        /** get str rep of gene in mGene1
         * @return string rep of gene in mGene1
         */
        inline std::string getGene1(){
            std::vector<std::string> trs;
            for(uint32_t i = 0; i < mGene1.size(); ++i){
                trs.push_back(mGene1[i].gene);
            }
            return util::join(trs, ";");
        }
        
        /** get str rep of gene in Gene2
         * @return string rep of gene in mGene2
         */
        inline std::string getGene2(){
            std::vector<std::string> trs;
            for(uint32_t i = 0; i < mGene2.size(); ++i){
                trs.push_back(mGene2[i].gene);
            }
            return util::join(trs, ";");
        }
        /** get str rep of transcripts of mGene1
         * @return string rep of transcript of mGene1
         */
        inline std::string getTrs1(){
            std::vector<std::string> trs;
            for(uint32_t i = 0; i < mGene1.size(); ++i){
                trs.push_back(mGene1[i].toStr());
            }
            return util::join(trs, ";");
        }
        
        /** get str rep of transcript of mGene2
         * @return string rep of transcript of mGene2
         */
        inline std::string getTrs2(){
            std::vector<std::string> trs;
            for(uint32_t i = 0; i < mGene2.size(); ++i){
                trs.push_back(mGene2[i].toStr());
            }
            return util::join(trs, ";");
        }
};

typedef std::vector<GeneInfo> GeneInfoList;
typedef std::vector<std::vector<BpRegion>> ContigBpRegions;
typedef std::vector<std::vector<SpanPoint>> ContigSpanPoints;

/** class to do coverage statistics of REF and ALT on contigs */
class Stats{
    public:
        Options* mOpt;                                     ///< pointer to Options
        std::vector<int32_t> mTotalAltCnts;                ///< total molecule(pair of read counted as one) supporting SV 
        std::vector<JunctionCount> mJctCnts;               ///< Single read spanning SV breakpoint stats
        std::vector<SpanningCount> mSpnCnts;               ///< Paired-end read spanning SV breakpoint stats

    public:
        /** Stats constructor */
        Stats(){}

        /** Stats constructor
         * @param opt pointer to Options object
         * @param n total SV number
         */
        Stats(Options* opt, int32_t n);

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
            for(uint32_t id = 0; id < st.mJctCnts.size(); ++id){
                os << "SV ID: " << id << "\n";
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
         * @param bpRegs SV breakpoint regions on each contig
         * @param spPts SV DP read mapping position on each contig
         * @param regInfo reference of RegItemCnt object
         * @param ctgCgrs sv events cgranges_t on each contig
         */
        void stat(const SVSet& svs, const ContigBpRegions& bpRegs, const ContigSpanPoints& spPts, const RegItemCnt& regInfo, cgranges_t* ctgCgr);

        /** merge coverage information of all contigs
         * @param sts reference of list of Stats
         * @param n number of svs in total
         * @param opt pointer to Options
         * @return merged stat info
         */
        static Stats* merge(const std::vector<Stats*>& sts, int32_t n, Options* opt);

        /** report BCF format report of all SVs
         * @param svs reference of SVSet
         */
        void reportSVBCF(const SVSet& svs);

        /** report TSV format report of all SVs
         * @param svs reference of SVSet
         * @param gl reference of GeneInfoList
         */
        void reportSVTSV(SVSet& svs, GeneInfoList& gl);

        /** report TSV format report of valid Fusion events
         * @param svs reference of SVSet
         * @param gl reference of GeneInfoList
         */
        void reportFusionTSV(const SVSet& svs, GeneInfoList& gl);
        
        /** convert an sv event to fusion record
         * @param fsr reference of FusionRecord
         * @param svr pointer of SVRecord
         * @param gi reference of GeneInfo
         * @param i which fusion will  be converted(-1 to all)
         */
        void toFuseRec(FusionRecord& fsr, const SVRecord* svr, GeneInfo& gi, int32_t i);

        /** mask fusion event status 
         * @param svs reference of SVSet
         * @param gl reference of GeneInfoList
         */
        void makeFuseRec(const SVSet& svs, GeneInfoList& gl);

        /** get alignment quality of sequence against an read
         * @param alnResult align result(sequence vertical, read horizontal)
         * @param qual read quality
         * @return average mapping quality of matched bases
         */
        uint32_t getAlignmentQual(Matrix2D<char>* alnResult, const uint8_t* qual);


        /** check whether an alignment of probe to read is valid
         * @param alnResult align result(sequence vertical, read horizontal)
         * @param bppos breakpoint position on read by ogirinam BAM aligner
         * @param seqlen length of read sequence
         * @param minoff minimum offset needed to spanning bppos
         * @return true if alnResult spanning bppos with at least minoff bps
         */
        bool validAlignment(Matrix2D<char>* alnResult, int32_t bppos,  int32_t seqlen, int32_t minoff = 0);

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
