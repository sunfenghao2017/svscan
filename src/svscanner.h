#ifndef SVSCANNER_H
#define SVSCANNER_H

#include "realnfilter.h"
#include "dpbamrecord.h"
#include "srbamrecord.h"
#include "annotator.h"
#include "svrecord.h"
#include "options.h"
#include "bamtotb.h"
#include <utility>
#include <vector>
#include <map>

/** class to scan sv */
class SVScanner{
    public:
        Options* mOpt;         ///< pointer to Options object
        RegionList mScanRegs; ///< valid regions to scan for SV
        SVSet mDPSVs;          ///< DP supported SV records
        SVSet mSRSVs;          ///< SR supported SVrecords
        ContigSRs mCtgSRs;     ///< SR supporting SV on each contig

    public:
        /** SVScanner constructor
         * @param opt pointer to Options
         */
        SVScanner(Options* opt){
            mOpt = opt;
            mScanRegs = opt->scanRegs;
        }

        /** SVScanner destructor */
        ~SVScanner(){}

        /** scan bam for DP and SR supporting SVs */
        void scanDPandSR();

        /** scan one contig for DP and SR supporting SVs
         * @param tid contig index to scanning DP and SR supporting SVs
         * @param jctMap to store SRs
         * @param dprSet to store DPs
         */
        void scanDPandSROne(int32_t tid, JunctionMap* jctMap, DPBamRecordSet* dprSet);

        /** test whether an read is valid to commpute SV by test its overlapping with creg and its mate overlapping with creg
         * @param b pointer to bam1_t
         * @param h pointer to bam_hdr_t
         * @return true if its valid
         */
        bool inValidReg(const bam1_t* b, const bam_hdr_t* h){
            if(!mOpt->overlapRegs) return true;
            if(mOpt->overlapRegs->overlap(h->target_name[b->core.tid], b->core.pos - b->core.l_qseq, b->core.pos + b->core.l_qseq)){
                return true;
            }
            if((b->core.flag & BAM_FMUNMAP) || b->core.mtid < 0){
                return false;
            }
            if(mOpt->overlapRegs->overlap(h->target_name[b->core.mtid], b->core.mpos - b->core.l_qseq - mOpt->libInfo->mReadLen, b->core.mpos + b->core.l_qseq + mOpt->libInfo->mReadLen)){
                return true;
            }
            uint8_t* sa = bam_aux_get(b, "SA");
            if(sa){
                std::string ss = bam_aux2Z(sa);
                std::vector<std::string> vsas;
                util::split(ss, vsas, ";");
                std::vector<std::string> vstr;
                for(uint32_t i = 0; i < vsas.size() - 1; ++i){
                    util::split(vsas[i], vstr, ",");
                    int32_t refpos = std::atoi(vstr[1].c_str()) - 1;
                    if(mOpt->overlapRegs->overlap(vstr[0].c_str(), refpos - b->core.l_qseq, refpos + b->core.l_qseq)){
                        return true;
                    }
                }
            }
            return false;
        }
};

#endif
