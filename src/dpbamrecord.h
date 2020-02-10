#ifndef DPBAMRECORD_H
#define DPBAMRECORD_H

#include <vector>
#include <cstdint>
#include <algorithm>
#include <htslib/sam.h>
#include "options.h"
#include "svrecord.h"

/** Class to store discordant pair of reads alignment record which supports SVs */
class DPBamRecord{
    public:
        int32_t mCurTid;    ///< current read(second read in coordinate sorted bam) bam record reference id b->core.tid
        int32_t mCurPos;    ///< current read(second read in coordinate sorted bam) bam record mapping starting pos b->core.pos
        int32_t mMateTid;   ///< mate read(first read in coordinate sorted bam) bam record reference id b->core.mtid
        int32_t mMatePos;   ///< mate read(first read in coordinate sorted bam) bam record mapping starting pos b->core.mpos
        int32_t mCurAlen;   ///< current read(second read in coordinate sorted bam) bam record consumed reference length bam_cigar2rlen(b)
        int32_t mMateAlen;  ///< mate read(first read in coordinate sorted bam) bam record consumed reference length bam_cigar2rlen(mateb)
        uint8_t mMapQual;   ///< b->core.qual, b is the record of read on large rcoordinate or larger chr in a pair
        int32_t mSVT;       ///< SV type this DPBamRecord supports
        int32_t mSVID = -1; ///< id of SV event this DPBamRecord supports;

    public:
        /** DPBamRecord constructor
         * @param b pointer to bam1_t struct(one read in a pair mapping at higher position or larger chromosome)
         * @param svt SV type this DPBamRecord supports
         */
        DPBamRecord(const bam1_t* b, int32_t svt){
            mCurTid = b->core.tid;
            mCurPos = b->core.pos;
            mMateTid = b->core.mtid;
            mMatePos = b->core.mpos;
            mCurAlen = bamutil::getSeqLen(b);
            std::vector<std::pair<int32_t, char>> pcigar;
            uint8_t* mctag = bam_aux_get(b, "MC");
            if(mctag){
                bamutil::parseCigar(bam_aux2Z(mctag), pcigar);
                mMateAlen = 0;
                for(auto& e: pcigar){
                    if(e.second == 'I' || e.second == 'X' || e.second == '=' || e.second == 'M' || e.second == 'S') mMateAlen += e.first;
                 }
            }else{
                mMateAlen = mCurAlen;
            }
            mMapQual = b->core.qual;
            mSVT = svt;
            mSVID = -1;
        }

        /** DPBamRecord destructor */
        ~DPBamRecord(){}

        /** operator to output an DPBamRecord to ostream
         * @param os reference of ostream object
         * @param dpr reference of DPBamRecord object
         * @return reference of ostream object
         */
        inline friend std::ostream& operator<<(std::ostream& os, const DPBamRecord& dpr){
            os << "===============================================================\n";
            os << "SVT: " << dpr.mSVT << "\n";
            os << "Current Read Tid: " << dpr.mCurTid << "\n";
            os << "Current Read Mapping Pos: " << dpr.mCurPos << "\n";
            os << "Current Read Reference Consumed: " << dpr.mCurAlen << "\n";
            os << "Mate Read Tid: " << dpr.mMateTid << "\n";
            os << "Mate Read Mapping Pos: " << dpr.mMatePos << "\n";
            os << "Mate Read Reference Consumed: " << dpr.mMateAlen << "\n";
            os << "Min Mapping Quality of Pair: " << dpr.mMapQual << "\n";
            os << "ID of Structural Variant this DP contributed to: " << dpr.mSVID << "\n";
            os << "===============================================================\n";
            return os;
        }

        /** operator to compare two DPBamRecord object\n
         * this compare rule will guarantee that DPBamRecord on same chr\n
         * will be sorted by their leftmost mapping position\n
         * DPBamRecord on different chr will be sorted by their mCurPos fields\n
         * @param other reference of another DPBamRecord object
         * @return true if this DPBamRecord object < other
         */
        inline bool operator<(const DPBamRecord& other) const {
            return (mCurTid < other.mCurTid) ||
                   (mCurTid == other.mCurTid && mMateTid < other.mMateTid) ||
                   (mCurTid == other.mCurTid && mMateTid == other.mMateTid && mCurPos < other.mCurPos) ||
                   (mCurTid == other.mCurTid && mMateTid == other.mMateTid && mCurPos == other.mCurPos && mMatePos < other.mMatePos);
        }

        /** return minimal coordinate of an DP\n
         * if DP is on same chr, then return the leftmost mapping position\n
         * else return mCurPos(the larger chr mapping pos)
         * @return minimal coordinate of DP
         */
        inline int32_t minCoord() const {
            if(mCurTid == mMateTid){
                return std::min(mCurPos, mMatePos);
            }else{
                return mCurPos;
            }
        }

        /** return maximum coordinate of an DP\n
         * if DP is on same chr, then return the right most mapping position\n
         * else return mCurPos(the smaller chr mapping pos)
         * @return maximum coordinate of DP
         */
        inline int32_t maxCoord() const {
            if(mCurTid == mMateTid){
                return std::max(mCurPos, mMatePos);
            }else{
                return mMatePos;
            }
        }
        
        /** get SV type of one pair of DP read supporting
         * @param b pointer to alignment record of one read of the DP pair
         * @param opt pointer to Options object
         */
        static int getSVType(const bam1_t* b, Options* opt);

        /** get SV type of one pair of DP read supporting
         * @param b pointer to alignment record of one read of the DP pair
         */
        static int getSVType(const bam1_t* b);

        /** initialize an clique(cluster of DPBamRecords which support the same SV event) wich an DPBamRecord
         * @param svStart to store SV starting position this DPBamRecord supporting
         * @param svEnd to store SV ending position this DPBamRecord supporting
         * @param svt type of SV event this DPBamRecord supprting [0-8]
         */
        void initClique(int32_t& svStart, int32_t& svEnd, int32_t svt);
        
        /** update an clique(cluster of DPBamRecords which support the same SV event) wich an DPBamRecord
         * @param svStart to store SV starting position this DPBamRecord supporting
         * @param svEnd to store SV ending position this DPBamRecord supporting
         * @param svt type of SV event this DPBamRecord supprting [0-8]
         */
        void updateClique(int32_t& svStart, int32_t& svEnd, int32_t svt);


};

/** class to store and analysis DPBamRecords */ 
class DPBamRecordSet{
    public:
        Options* mOpt;///< pointer to Options object
        std::vector<std::vector<DPBamRecord>> mDPs;///< DPBamRecord supporting various types of SVs
    
    public:
        /** DPBamRecordSet constructor
         * @param opt pointer to Options object
         */
        DPBamRecordSet(Options* opt){
            mOpt = opt;
            mDPs.resize(9);
        }

        /** DPBamRecordSet destructor */
        ~DPBamRecordSet(){
            mOpt = NULL;
        }

    public:
        /** merge a list of DPBamRecordSet into one
         * @param dps list of DPBamRecordSet
         * @param opt pointer to Options
         * @return merged DPBamRecordSet
         */
        static inline DPBamRecordSet* merge(const std::vector<DPBamRecordSet*> dps, Options* opt){
            DPBamRecordSet* ret = new DPBamRecordSet(opt);
            for(auto& e: dps){
                for(int32_t i = 0; i < 9; ++i){
                    std::copy(e->mDPs[i].begin(), e->mDPs[i].end(), std::back_inserter(ret->mDPs[i]));
                }
            }
            return ret;
        }

        /** operator to output an DPBamRecordSet to ostream
         * @param os reference of ostream object
         * @param dps reference of DPBamRecordSet object
         * @return reference of ostream object
         */
        inline friend std::ostream& operator<<(std::ostream& os, const DPBamRecordSet& dps){
            for(uint32_t i = 0; i < dps.mDPs.size(); ++i){
                os << "SVT: " << i << ":" << dps.mDPs[i].size() << "\n";
                for(uint32_t j = 0; j < dps.mDPs[i].size(); ++j){
                    os << dps.mDPs[i][j];
                }
                os << "\n";
            }
            return os;
        }
        
        /** operator to output an DPBamRecordSet to ostream
         * @param os reference of ostream object
         * @param dps reference of DPBamRecordSet object
         * @return reference of ostream object
         */
        inline friend std::ostream& operator<<(std::ostream& os, const DPBamRecordSet* dps){
            for(uint32_t i = 0; i < dps->mDPs.size(); ++i){
                os << "SVT: " << i << ":" << dps->mDPs[i].size() << "\n";
                for(uint32_t j = 0; j < dps->mDPs[i].size(); ++j){
                    os << dps->mDPs[i][j];
                }
                os << "\n";
            }
            return os;
        }

        /** inert an DPBamRecord
         * @param b pointer to bam1_t* struct which is the latter read in coordinate sorted bam which is DP
         * @param svt SV type this bam record supports
         */
        inline void insertDP(const bam1_t* b, int32_t svt){
            mDPs[svt].push_back(DPBamRecord(b, svt));
        }

        /** cluster sorted DPBamRecords of one type SV and find all supporting SV of this type\n
         * step1: cluster DPBamRecord into different component, DPBamRecord with same chr1 and starting and ending\n
         *        mapping positions are in reasonable range(only considering nearing two) consists a component
         * step2: search each component for an clique, use edge with smallest weight as seed,\n
         *        this edge have the smallest leftmost and rightmost coordinates offset between src and dest\n
         *        initialize SV starting/ending wiggling(maximum mapping pos offset beetween an clique) in a\n
         *        SV type specific way, then expand this clique as big as possible\n
         * @param dps a list of DPBamRecord which supporting one kind of SV type
         * @param svs SVSet to store SV supporting by DPs
         * @param svt only cluster this type of SV supporting DPBamRecords
         */
        void cluster(std::vector<DPBamRecord>& dps, SVSet* svs, int32_t svt);
        
        /** cluster all kinds SV supporting DPBamRecord
         * @param svs SVSet which support various kind of SVs
         */
        void cluster(SVSet* svs);
        
        /** a subroutine used to search all possible clique supporting an type of SV
         * @param clique discordant pairs clustered from DPBamRecords
         * @param dps reference of DPBamRecordSet list which supporting one kind of SV
         * @param svs SVSet to store SV supporting by DPs
         * @param svt SV type analyzed
         */
        void searchCliques(std::set<int32_t>& clique, std::vector<DPBamRecord>& dps, SVSet* svs, int32_t svt);
       
        /** check whether SV size is valid 
         * @param svStart SV starting position
         * @param svEnd SV ending position
         * @param svt SV type
         * @return true if valid
         */ 
        inline bool validSVSize(int32_t svStart, int32_t svEnd, int32_t svt){
            if(svt == 0 || svt == 1) return svEnd - svStart >= mOpt->filterOpt->mMinInversionRpt;
            if(svt == 2) return svEnd - svStart >= mOpt->filterOpt->mMinDeletionRpt;
            if(svt == 3) return svEnd - svStart >= mOpt->filterOpt->mMinDupRpt;
            return true;
        }
};
#endif
