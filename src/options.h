#ifndef OPTIONS_H
#define OPTIONS_H

#include <set>
#include <mutex>
#include <cstdio>
#include <vector>
#include <string>
#include <cstdint>
#include <sstream>
#include "realnfilter.h"
#include <htslib/sam.h>
#include "fusionopt.h"
#include "ThreadPool.h"
#include "software.h"
#include "statutil.h"
#include "util.h"
#include "bed.h"

#define DEBUG_FCALL 0x1  ///< calling debug mask
#define DEBUG_FANNC 0x2  ///< coverage anno debug mask
#define DEBUG_FANNG 0x4  ///< gene anno debug mask
#define DEBUG_FOUTF 0x8  ///< output debug mask
#define DEBUG_FREAN 0x10 ///< realign debug mask

typedef uint32_t DEBUG_TYPE; ///< debug type

/** class to store library information */
struct LibraryInfo{
    int32_t mMaxSample = 1000000; ///< maximum number of reads used to generate library information
    int32_t mReadLen = 0;         ///< read length
    int32_t mMedian = 0;          ///< isize median
    int32_t mMad = 0;             ///< isize mad (median of abs(isize - median)
    int32_t mMinNormalISize = 0;  ///< normal minimum isize max(0, median - 5 * mad)
    int32_t mMinISizeCutoff = 0;  ///< customized isize minimum cutoff
    int32_t mMaxNormalISize = 0;  ///< normal maximum isize (median + 5 * mad)
    int32_t mMaxISizeCutoff = 0;  ///< max(0, median + cutoff * mad, 2 * readlen, 500)
    int32_t mAbnormalPairs = 0;   ///< abnormal discordant pairs count
    bool mIsHaploTagged = false;  ///< bam has HP tag if true
    int32_t mContigNum = 0;       ///< maximum contig number in library
    int32_t mVarisize = 0;        ///< std::max(mMaxNormalISize, mReadLen), mMaxDPVarSize equqls this value
    
    /** LibraryInfo constructor */
    LibraryInfo(){}
    
    /** LibraryInfo destructor */
    ~LibraryInfo(){}

    /** convert libInfo to string
     * @return string representation of LibraryInfo
     */
    std::string toStr(){
        std::stringstream ss;
        ss << "Read Length: " << mReadLen << "\n";
        ss << "ISize Median: " << mMedian << "\n";
        ss << "ISize MAD: " << mMad << "\n";
        ss << "Minimum Normal ISize: " << mMinNormalISize << "\n";
        ss << "Minimum ISize cutoff: " << mMinISizeCutoff << "\n";
        ss << "Maximum Normal ISize: " << mMaxNormalISize << "\n";
        ss << "Maximum ISize cutoff: " << mMaxISizeCutoff << "\n";
        ss << "Abnormal DP read: " << mAbnormalPairs << "\n";
        ss << "Library Haplotype Tagged: " << std::boolalpha << mIsHaploTagged << "\n";
        ss << "Contig/Chrosome Number: " << mContigNum << "\n";
        ss << "Maximum Varsize : " << mVarisize;
        return ss.str();
    }

    /** operator to output LibraryInfo to ostream
     * @param os reference of ostream
     * @param libinfo referencce of LibraryInfo
     * @return reference of ostream
     */
    friend std::ostream& operator<<(std::ostream& os, const LibraryInfo& libinfo){
        os << "Read Length: " << libinfo.mReadLen << "\n";
        os << "ISize Median: " << libinfo.mMedian << "\n";
        os << "ISize MAD: " << libinfo.mMad << "\n";
        os << "Minimum Normal ISize: " << libinfo.mMinNormalISize << "\n";
        os << "Minimum ISize cutoff: " << libinfo.mMinISizeCutoff << "\n";
        os << "Maximum Normal ISize: " << libinfo.mMaxNormalISize << "\n";
        os << "Maximum ISize cutoff: " << libinfo.mMaxISizeCutoff << "\n";
        os << "Abnormal DP read: " << libinfo.mAbnormalPairs << "\n";
        os << "Library Haplotype Tagged: " << std::boolalpha << libinfo.mIsHaploTagged << "\n";
        os << "Contig/Chrosome Number: " << libinfo.mContigNum << "\n";
        os << "Maximum Varsize : " << libinfo.mVarisize << "\n";
        return os;
    }
    
    /** operator to output LibraryInfo to ostream
     * @param os reference of ostream
     * @param libinfo pointer to LibraryInfo
     * @return reference of ostream
     */
    friend std::ostream& operator<<(std::ostream& os, LibraryInfo* libinfo){
        os << *libinfo << std::endl;
        return os;
    }
};

/** class to store various filter options to SVs or SV supporting SR/DP reads */
struct SVFilter{
    int32_t mMinRefSep = 50;           ///< minimal reference seperation needed for an split alignment used to compute SV
    int32_t mMaxReadSep = 10;          ///< maximal read split alignment position(both from 5') allowed to be used to compute SV
    double mFlankQuality = 0.95;       ///< flank identity ratio...
    int32_t mMinFlankSize = 10;        ///< minimal flank length needed for consensus split read length on each side of breakpoint
    int32_t minMapQual = 1;            ///< minimal paired-end(PE) mapping quality
    int32_t minClipLen = 20;           ///< minimal clipping length used to compute SV
    int32_t mMinTraQual = 1;           ///< minimal PE mapping quality for translocation
    int32_t mMinDupSize = 100;         ///< minimal duplication size needed for an DP record used to compute SV
    uint32_t mMinGenoQual = 1;         ///< minimal mapping quality for genotyping
    uint32_t mGraphPruning = 100000;   ///< PE graph pruning cutoff
    int32_t mMaxReadPerSV = 100000;    ///< maximum valid split-reads sampled to analysis one SV event
    int32_t mMinGapOfCSSVRef = 15;     ///< minimal middle gap length needed for an SR consensus against SV ref seq alignment
    int32_t mMaxCoordDevOfCSSVRef = 5; ///< maximum gap start/end length allowed for the non-gapped partner of an valid SR consensus ~ SV ref seq alignment
    int32_t mMinInversionRpt = 100;    ///< minimum inversion size to report
    int32_t mMinDeletionRpt = 300;     ///< minimum deletion size to report
    int32_t mMinDupRpt = 100;          ///< minimum duplication size to report
    int32_t mMinInsRpt = 15;           ///< minimum insertion size to repor
    uint32_t mMinSeedSR = 3;           ///< minimum seed split reads used to compute SV
    uint32_t mMinSeedDP = 3;           ///< minimum discordant pair of reads used to compute SV
    float mMinDelRatio = 0.8;          ///< minimum deletion ratio of an exon to report 
    float mMaxFPIns = 0.5;             ///< maximum ratio of reads supporting both INS and other type SVs allowed for an valid insertion
    float mMinSRResScore = 1;          ///< minimal alignment score for SR rescued

    /** SVFilter constructor */
    SVFilter(){}

    /** SVFilter destructor */
    ~SVFilter(){}
};

/** class to store threshold of high quality SV */
struct PassOptions{
    int32_t mIntraChrSVMinDPCnt = 2;   ///< intra-chromosome SV minimum DP supports
    int32_t mIntraChrSVMinDPQual = 15; ///< intra-chromosome SV minimum DP quality
    int32_t mIntraChrSVMinSRCnt = 2;   ///< intra-chromosome SV minimum SR supports
    int32_t mIntraChrSVMinSRQual = 15; ///< intra-chromosome SV minimum DP quality
    int32_t mInterChrSVMinDPCnt = 2;   ///< inter-chromosome SV minimum DP supports
    int32_t mInterChrSVMinDPQual = 15; ///< inter-chromosome SV minimum DP quality
    int32_t mInterChrSVMinSRCnt = 2;   ///< inter-chromosome SV minimum SR supports
    int32_t mInterChrSVMinSRQual = 15; ///< inter-chromosome SV minimum DP quality
    int32_t mMinGenoQual = 15;         ///< minimum genotype quality

    /** PassOptions constructor */
    PassOptions(){}

    /** PassOptions destructor */
    ~PassOptions(){}
};

/** class to store MSA argument */
struct MSAOpt{
    int32_t mMinCovForCS = 3;          ///< minimum coverage needed for a position in msa result to be included in consensus sequence
    float mMinBaseRateForCS = 0.5;     ///< minimum base ratio needed for a position in msa result to be included in consensus sequence
    bool mAlignHorzEndGapFree = false; ///< use horizontal end gap penalty free strategy when get consensus sequence of SR
    bool mALignVertEndGapFree = false; ///< use vertical end gap penalty free strategy when get consensus sequence of SR

    /** MSAOpt consturcor */
    MSAOpt(){}

    /** MSAOpt destructor */
    ~MSAOpt(){}
};

/** type to store regions */
typedef std::vector<std::set<std::pair<int32_t, int32_t>>> RegionList;

/** class to store various options */
class Options{
    public:
        std::string bamfile;          ///< bam file used to analysis currently
        std::string alnref;           ///< reference used to align reads in bam
        std::string genome;           ///< genome reference of the species(rna mode)
        std::string annodb;           ///< annotation feature database file
        std::string gannodb;          ///< genome annotation database file(rna mode)
        std::string reg;              ///< file to store regions to scanning bam in
        std::string creg;             ///< file to store regions that SV event must overlap
        std::string bcfOut;           ///< output SV bcf result file
        std::string tsvOut;           ///< output SV tab seperated values file
        std::string bamout;           ///< output SV supporting bam record file
        std::string bam2tb;           ///< fusion event supporting reads excel file
        samFile* fbamout;             ///< file pointer of sv bam output
        int32_t madCutoff;            ///< insert size cutoff, median+s*MAD (deletions only)
        int32_t batchsvn;             ///< batch sv events to annotate coverage 
        RegionList scanRegs;          ///< regions to scan bam
        BedRegs* overlapRegs;         ///< regions sv must overlap
        std::vector<int32_t> svtypes; ///< sv types to discovery(for commandline argument parsing)
        std::set<int32_t> SVTSet;     ///< predefined sv types to compute [INV, DEL, DUP, INS, BND]
        int32_t nthread;              ///< threads used to process REF/ALT read/pair assignment
        ThreadPool::ThreadPool* pool; ///< thread pool used to do parallel work
        std::mutex logMtx;            ///< mutex locked to output log information
        std::mutex outMtx;            ///< mutex locked to output sv supporting bam
        int32_t contigNum;            ///< max contig numbers in library bam
        std::set<int32_t> svRefID;    ///< SV occuring reference id
        LibraryInfo* libInfo;         ///< library information for the currently analyzed bam
        SVFilter* filterOpt;          ///< filter options
        PassOptions* passOpt;         ///< high quality SV threshold
        FusionOptions* fuseOpt;       ///< fusion report options
        MSAOpt* msaOpt;               ///< MSA options
        Software* softEnv;            ///< softare environment option
        bam_hdr_t* bamheader;         ///< bam header
        RealnFilter* realnf;          ///< realign filter
        bool rnamode;                 ///< find rna structural variants
        bool writebcf;                ///< write bcf if true
        DEBUG_TYPE debug;             ///< debug mode

    public:
        /** Options constructor */
        Options();

        /** Options destructor */
        ~Options();

        /** validate some arguments passed */
        void validate();

        /** update options after command arguments parsed
         * @param argc number of arguments feed to main
         * @param argv array of arguments feed to main
         */
        void update(int argc, char** argv);

        /** update library informations
         * @param bam bam file path
         */
        LibraryInfo* getLibInfo(const std::string& bam);

        /** create valid regions by exclude invalid regions */
        void getScanRegs();

        /** get regions each sv event must overlap */
        void getCregs();

        /** write empty result file in case input bam is empty */
        void writeEmptFile();
};

#endif
