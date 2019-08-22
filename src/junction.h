#ifndef JUNCTION_H
#define JUNCTION_H

#include <map>
#include <string>
#include <cstdint>
#include <algorithm>
#include "svutil.h"
#include "options.h"
#include <htslib/sam.h>

/** class to store junction alignment record */
class Junction{
    public:
        bool mForward;     ///< junction read is from forward strand if true (eg. !b->core.BAM_FREVERSE)
        bool mSCleft;      ///< softclip is at leading left part of alignment if true
        int32_t mSCLen;    ///< softclip length of this junction alignment record
        int32_t mRefidx;   ///< junction record alignment reference tid (b->core.tid)
        int32_t mRstart;   ///< junction record alignment starting position on reference(b->core.pos)
        int32_t mRefpos;   ///< b->core.pos + reference length consumed before junction point
        int32_t mSeqpos;   ///< sequence length consumed before junction point(count from read 5'->3')
        int32_t mSeqmatch; ///< sequence length consumed before junction point(count from read alignment direction)
    public:
        /** Junction object constructor 
         * @param forward junction read is from forward strand if true
         * @param scleft softclip is at leading left part of alignment if true
         * @param sclen softclip length
         * @param refidx junction record alignment reference tid (b->core.tid)
         * @param rstart junction record alignment starting position on reference(b->core.pos)
         * @param refpos b->core.pos + reference length consumed before junction point
         * @param seqpos sequence length consumed before junction point(count from read 5'->3')
         * @param seqmatch sequence length consumed before junction point(count from read alignment direction)
         */
        Junction(bool forward, bool scleft, int32_t sclen, int32_t refidx, int32_t rstart, int32_t refpos, int32_t seqpos, int32_t seqmatch){
            mForward = forward;
            mSCleft = scleft;
            mSCLen = sclen;
            mRefidx = refidx;
            mRstart = rstart;
            mRefpos = refpos;
            mSeqpos = seqpos;
            mSeqmatch = seqmatch;
        }

        /** operator to output an Junction object to ostream
         * @param os reference of ostream object
         * @param jct reference of Junction object
         * @return reference of ostream object
         */
        inline friend std::ostream& operator<<(std::ostream& os, const Junction& jct){
            os << "==========================================\n";
            os << std::boolalpha << "From Forward Strand: " << jct.mForward << "\n";
            os << std::boolalpha << "Leading Soft Clip: " << jct.mSCleft << "\n";
            os << "Soft clip length: " << jct.mSCLen << "\n";
            os << "Reference ID: " << jct.mRefidx << "\n";
            os << "Reference Mapping Pos: " << jct.mRstart << "\n";
            os << "Clip position on Ref: " << jct.mRefpos << "\n";
            os << "Clip position on Read: " << jct.mSeqpos << "\n";
            os << "Length eat before/after Clip: " << jct.mSeqmatch << "\n";
            os << "==========================================\n";
            return os;
        }

        /** operator to compare two Junction object
         * @return true if this Junction object is less than other
         */
        inline bool operator<(const Junction& other) const {
            return mSeqpos < other.mSeqpos ||
               (mSeqpos == other.mSeqpos && mRefidx < other.mRefidx) ||
               (mSeqpos == other.mSeqpos && mRefidx == other.mRefidx && mRefpos < other.mRefpos) ||
               (mSeqpos == other.mSeqpos && mRefidx == other.mRefidx && mRefpos == other.mRefpos && mSCleft < other.mSCleft);
        }
};

/** Class to store Junction reads */
class JunctionMap{
    public:
        Options* mOpt;                                          ///< pointer to Options
        std::map<size_t, std::vector<Junction>> mJunctionReads; ///< key:hash value of read name, value: junction read part list
        bool mSorted;                                           ///< Junction records in mJunctionReads are all sorted if true

    public:
        /** JunctionMap constructor
         * @param opt pointer to Options
         */
        JunctionMap(Options* opt){
            mOpt = opt;
            mSorted = false;
        }

        /** JunctionMap destructor */
        ~JunctionMap(){}
    public:

        /** insert an read to JunctionMap if it is junction read
         * @param b pointer to bam1_t struct
         * @param h pointer to bam_hdr_t struct
         * @return true if b is an junction read and inserted successfuly
         */
        bool insertJunction(const bam1_t* b, bam_hdr_t* h);

        /** sort all Junction records in mJunctionReads */
        void sortJunctions();
        
        /** operator to output an JunctionMap object to ostream
         * @param os reference of ostream object
         * @param jctMap reference of JunctionMap object
         * @return reference of ostream object
         */
        inline friend std::ostream& operator<<(std::ostream& os, const JunctionMap& jctMap){
            for(auto& e: jctMap.mJunctionReads){
                os << "Read Hash: " << e.first << "\n";
                int i = 1;
                for(auto& f: e.second) os << "===== " << i++ << " =====\n" << f;
                os << "\n";
            }
            return os;
        }

        /** merge a list of JunctionMap into one
         * @param jct a list of JunctionMap
         * @param opt pointer to Options
         * @return merged JunctionMap
         */
        static inline JunctionMap* merge(const std::vector<JunctionMap*>& jct, Options* opt){
            JunctionMap* ret = new JunctionMap(opt);
            for(auto& e: jct){
                std::copy(e->mJunctionReads.begin(), e->mJunctionReads.end(), std::inserter(ret->mJunctionReads, ret->mJunctionReads.end()));
            }
            return ret;
        }
};

#endif
