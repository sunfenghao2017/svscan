#ifndef FUSION_REPORTER_H
#define FUSION_REPORTER_H

#include "fusionopt.h"
#include "software.h"
#include "CLI.hpp"

struct FusionRecord{
    std::string fusegene;
    std::string fusepattern;
    int32_t fusionreads;
    int32_t totalreads;
    float fuserate;
    std::string gene1;
    std::string chr1;
    int32_t junctionposition1;
    std::string strand1;
    std::string transcript1;
    std::string gene2;
    std::string chr2;
    int32_t junctionposition2;
    std::string strand2;
    std::string transcript2;
    std::string fusionsequence;
    std::string indb;
    std::string svt;
    std::string svsize;
    std::string srcount;
    std::string dpcount;
    std::string srrescued;
    std::string dprescued;
    std::string srrefcount;
    std::string dprefcount;
    std::string insbp;
    std::string insseq;
    std::string svid;
    std::string svint;
    bool report;

    FusionRecord(){}

    ~FusionRecord(){}

    inline friend std::ostream& operator<<(std::ostream& os, const FusionRecord& fsr){
        os << fsr.fusegene << "\t" << fsr.fusepattern << "\t";
        os << fsr.fusionreads << "\t" << fsr.totalreads << "\t" << fsr.fuserate << "\t";
        os << fsr.gene1 << "\t" << fsr.chr1 << "\t" << fsr.junctionposition1 << "\t" << fsr.strand1 << "\t" << fsr.transcript1 << "\t";
        os << fsr.gene2 << "\t" << fsr.chr2 << "\t" << fsr.junctionposition2 << "\t" << fsr.strand2 << "\t" << fsr.transcript2 << "\t";
        os << fsr.fusionsequence << "\t" << fsr.indb << "\t" << fsr.svt << "\t" << fsr.svsize << "\t";
        os << fsr.srcount << "\t" << fsr.dpcount << "\t" << fsr.srrescued << "\t" << fsr.dprescued << "\t";
        os << fsr.srrefcount << "\t" << fsr.dprefcount << "\t";
        os << fsr.insbp << "\t" << fsr.insseq << "\t" << fsr.svid << "\t" << fsr.svint << "\n";
       return os;
    } 
};

typedef std::vector<FusionRecord> FusionList;

struct FusionReporter{
    FusionOptions* fuseOpt; ///< fusion report options
    Software* softEnv;      ///< software environments
    FusionList fuseList;    ///< fusion list

    /** FusionReporter constructor */
    FusionReporter();

    /** FusionReporter destructor */
    ~FusionReporter();

    /** update options after command arguments parsed
    * @param argc number of arguments feed to main
    * @param argv array of arguments feed to main
    */
    void update(int argc, char** argv);
   
    /** parse sv.tsv into FusionList */
    void sv2fs(); 

    /** do fusion report job */
    void report();
};

#endif
