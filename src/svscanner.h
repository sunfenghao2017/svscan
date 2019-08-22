#ifndef SVSCANNER_H
#define SVSCANNER_H

#include "annotator.h"
#include "dpbamrecord.h"
#include "srbamrecord.h"
#include "svrecord.h"
#include "options.h"
#include <utility>
#include <vector>
#include <map>

/** class to scan sv */
class SVScanner{
    public:
        Options* mOpt;         ///< pointer to Options object
        RegionList mValidRegs; ///< valid regions to scan for SV
        SVSet mDPSVs;          ///< DP supported SV records
        SVSet mSRSVs;          ///< SR supported SVrecords
        ContigSRs mCtgSRs;     ///< SR supporting SV on each contig

    public:
        /** SVScanner constructor
         * @param opt pointer to Options
         */
        SVScanner(Options* opt){
            mOpt = opt;
            mValidRegs = opt->validRegions;
        }

        /** SVScanner destructor */
        ~SVScanner(){}

        /** scan bam for DP and SR supporting SVs */
        void scanDPandSR();

        /** scan one conrig for DP and SR supporting SVs
         * @param tid contig index to scanning DP and SR supporting SVs
         * @param jctMap to store SRs
         * @param dprSet to store DPs
         */
         void scanDPandSROne(int32_t tid, JunctionMap* jctMap, DPBamRecordSet* dprSet);
};

#endif
