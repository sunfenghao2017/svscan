#ifndef FUSION_REPORTER_H
#define FUSION_REPORTER_H

#include "fusionopt.h"
#include "software.h"
#include "fuserec.h"
#include "svutil.h"
#include "trsrec.h"

struct FusionReporter{
    FusionOptions* fuseOpt; ///< fusion report options
    Software* softEnv;      ///< software environments
    FusionRecordList fuseList;    ///< fusion list
    bool rnamode;           ///< fusion event is from rnaseq

    /** FusionReporter constructor */
    FusionReporter();

    /** FusionReporter destructor */
    ~FusionReporter();

    /** update options after command arguments parsed
    * @param argc number of arguments feed to main
    * @param argv array of arguments feed to main
    */
    void update(int argc, char** argv);

    /** parse string to TrsRecList
     * @param trsl TrsRecList to store transcripts
     * @param trStr bp1Gene or bp2Gene field in sv.tsv
     */
    void str2trsl(TrsRecList& trsl, const std::string& trStr);

    /** parse string to FuseGeneList
     * @param fsgs FusionGene store
     * @param fsstr fuseGene field in sv.tsv
     * @param bp1trs TrsRecList of bp1Gene in sv.tsv
     * @param bp2trs TrsRecList of bp2Gene in sv.tsv
     */
    void str2fsgs(FuseGeneList& fsgl, const std::string& fsStr, const TrsRecList& bp1trs, const TrsRecList& bp2trs);
   
    /** parse sv.tsv into FusionRecordList */
    void sv2fs(); 

    /** do fusion report job */
    void report();
};

#endif
