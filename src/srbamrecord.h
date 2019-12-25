#ifndef SRBAMRECORD_H
#define SRBAMRECORD_H

#include <vector>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <htslib/faidx.h>
#include "msa.h"
#include "bamutil.h"
#include "options.h"
#include "junction.h"
#include "svrecord.h"

/** class to store split read alignment record */
class SRBamRecord{
    public:
        int32_t mChr1;     ///< reference id part1 of read mapped
        int32_t mPos1;     ///< break point position of part1 read on reference
        int32_t mChr2;     ///< reference id part2 of read mapped
        int32_t mPos2;     ///< break point position of part2 read on reference
        int32_t mRstart;   ///< starting mapping position on reference of one part read which is not -1
        int32_t mInslen;   ///< insert size of two part of one read contributed, this might caused by insertion or sequence error
        int32_t mSVID;     ///< default -1, if allocated to an StructuralVariant, it is the index at which to store a StructuralVariant in vector
        size_t mID;        ///< hash value of the read name
        bool mRead1;       ///< from read1 if true

    public:
        /** construct SRBamRecord object
         * @param chr1 reference id part1 of read mapped
         * @param pos1 break point position of part1 read on reference
         * @param chr2 reference id part2 of read mapped
         * @param pos2 break point position of part2 read on reference
         * @param rstart starting mapping position on reference of one part read which is not -1
         * @param inslen insert size of two part of one read contributed
         * @param id hash value of the read name
         * @param rd1 read1 if true
         */
        SRBamRecord(int32_t chr1, int32_t pos1, int32_t chr2, int32_t pos2, int32_t rstart, int32_t inslen, size_t id, bool rd1){
            mChr1 = chr1;
            mPos1 = pos1;
            mChr2 = chr2;
            mPos2 = pos2;
            mRstart = rstart;
            mInslen = inslen;
            mSVID = -1;
            mID = id;
            mRead1 = rd1;
        }

        /** SRBamRecord destructor */
        ~SRBamRecord(){}
        
        /** adjust orientation of SR read sequence to restore natural sequence in sample, only needed in 5to5 and 3to3 catenation
         * @param seq SR bam record seq
         * @param bpPoint if(SR mapped on little chr in translocation || mapped on higher coordinate in inversion) bpPoint = true;
         * @param svt SV type
         */
        inline static void adjustOrientation(std::string& seq, bool bpPoint, int32_t svt){
            if((svt == 5 && bpPoint) || (svt == 6 && !bpPoint) ||
               (svt == 0 && bpPoint) || (svt == 1 && !bpPoint)){
                util::reverseComplement(seq);
            }
        }

        /** operator to output an SRBamRecord object to ostream
         * @param os reference of ostream object
         * @param sr reference of SRBamRecord object
         * @return reference of ostream
         */
        inline friend std::ostream& operator<<(std::ostream& os, const SRBamRecord& sr){
            os << "===============================================================\n";
            os << "Split read is from read1 of pair: " << std::boolalpha << sr.mRead1 << "\n";
            os << "Part1 of Split Read Reference ID: " << sr.mChr1 << "\n";
            os << "Part1 of Split Read Breakpoint Position on Reference: " << sr.mPos1 << "\n";
            os << "Part2 of Split Read Reference ID: " << sr.mChr2 << "\n";
            os << "Part2 of SPlit Read Breakpoint Position on Reference: " << sr.mPos2 << "\n";
            os << "Insert size contribured by Part1 and Part2 of Split Read: " << sr.mInslen << "\n";
            os << "ID of Structural Variant this Split Read contributed to: " << sr.mSVID << "\n";
            os << "Hash value of Split Read name: " << sr.mID << "\n";
            os << "===============================================================\n";
            return os;
        }

        /** operator to compare two SRBamRecord object
         * @param other reference of SRBamRecord
         * @return true if this SRBamRecord is less than other
         */
        inline bool operator<(const SRBamRecord& other) const {
            return (mChr1 < other.mChr1) ||
                   (mChr1 == other.mChr1 && mChr2 < other.mChr2) || 
                   (mChr1 == other.mChr1 && mChr2 == other.mChr2 && mPos1 < other.mPos1) ||
                   (mChr1 == other.mChr1 && mChr2 == other.mChr2 && mPos1 == other.mPos1 && mPos2 < other.mPos2);
        }
};

/** class to store SR information on each contig */
typedef std::vector<std::map<std::pair<int32_t, size_t>, int32_t>> ContigSRs; ///<[contig]<<mRstart, mID>, mSVID>

/** class to store SRBamRecord supporting various SVs */
class SRBamRecordSet{
    public:
        Options* mOpt;                                        ///< pointer to Options object
        std::vector<std::vector<SRBamRecord>> mSRs;           ///< vector to store SRBamRecords according to the SV they support
        ContigSRs mSRMapPos;                                  ///< SR mapping starting position on each contig with SV Type defined
        bool mSorted = false;                                 ///< all SRBamRecords have been sorted if true
        std::vector<std::multiset<std::string>> mTraSeqStore; ///< translocation SR read sequence
        std::vector<std::multiset<std::string>> mTriSeqStore; ///< translocation insertion sequence nearby bp
        std::vector<std::vector<uint8_t>> mTraQualStore;      ///< translocation SR read mapping quality
        std::mutex mLock;                                     ///< lock used to update translocation info

    public:
        /** SRBamRecordSet constructor 
         * @param opt pointer to Options
         * @param jctMap pointer to JunctionMap
         */
        SRBamRecordSet(Options* opt, JunctionMap* jctMap = NULL){
            mOpt = opt;
            mSRs.resize(9);
            mSRMapPos.resize(mOpt->contigNum);
            if(jctMap) classifyJunctions(jctMap);
        }

        /** SRBamRecordSet destructor */
        ~SRBamRecordSet(){}

        /** sort SRBamRecords in mSRs */
        void sortSRs(){
            for(auto& svt : mOpt->SVTSet){
                std::sort(mSRs[svt].begin(), mSRs[svt].end());
            }
            mSorted = true;
        }
        
        /** operator to output an SRBamRecordSet object to ostream
         * @param os reference of ostream object
         * @param srs reference of SRBamRecordSet object
         * @return reference of ostream
         */
        inline friend std::ostream& operator<<(std::ostream& os, const SRBamRecordSet& srs){
            for(uint32_t i = 0; i < srs.mSRs.size(); ++i){
                os << "SVT: " << i << "\n";
                for(uint32_t j = 0; j < srs.mSRs[i].size(); ++j){
                    os << "===== " << j << " =====\n";
                    os << srs.mSRs[i][j];
                }
                os << "\n";
            }
            return os;
        }

        /** operator to output an SRBamRecordSet object to ostream
         * @param os reference of ostream object
         * @param srs pointer to SRBamRecordSet object
         * @return reference of ostream
         */
        inline friend std::ostream& operator<<(std::ostream& os, const SRBamRecordSet* srs){
            for(uint32_t i = 0; i < srs->mSRs.size(); ++i){
                os << "SVT: " << i << "\n";
                for(uint32_t j = 0; j < srs->mSRs[i].size(); ++j){
                    os << "===== " << j << " =====\n";
                    os << srs->mSRs[i][j];
                }
                os << "\n";
            }
            return os;
        }

        /** class all Junction reads record in JunctionMap into SRBamRecord vector according to SV type they support
         * @param jctMap pointer to JunctionMap
         */
        void classifyJunctions(JunctionMap* jctMap);

        /** cluster SRBamRecord of one type SV and find all supporting SV of this type\n
         * @param srs reference of SRBamRecords which supporting SV type svt
         * @param svs SVSet used to store SV found
         * @param svt SV type to find in srs, range [0-8]
         */
        void cluster(std::vector<SRBamRecord>& srs, SVSet& svs, int32_t svt);
        
        /** cluster all SRBamRecord in SRBamRecordSet into their seperate supporting SVs */
        void cluster(SVSet& svs);

        /** get sv candidates of a cluster
         * @param clique sr id of a cluster
         * @param srs reference of SRBamRecords which supporting SV type svt
         * @param svs SVSet used to store SV found
         * @param svt SV type to find in srs, range [0-8]
         */ 
        void searchCliques(std::set<int32_t>& clique, std::vector<SRBamRecord>& srs, SVSet& svs, int32_t svt);

        /** assembly reads of SR supporting each SV by MSA to get an consensus representation of SRs,\n
         * split align the consensus sequence against the constructed reference sequence to refine the breakpoint position
         * @param svs reference of SVSet
         */
        void assembleSplitReads(SVSet& svs);

        void assembleOneContig(SVSet& svs, int32_t refIdx);

        void assembleCrossChr(SVSet& svs, AlignConfig* alnCfg, bam_hdr_t* hdr, const std::vector<int32_t>& crsidx, int32_t begIdx, int32_t endIdx);
};

#endif
