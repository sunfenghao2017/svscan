#ifndef FUSEREC_H
#define FUSEREC_H

#include <string>
#include <cstdint>
#include <iostream>
#include "svutil.h"

/** class to store an fusion record */
struct FusionRecord{
    std::string fusegene;       ///< fusion gene hgene->tgene
    std::string fusepattern;    ///< fusion gene strand e.g., +->-, -->-, +->+
    int32_t fusionreads;        ///< # of reads supporting fusion event
    int32_t totalreads;         ///< # of reads supporting fusion event and refrence type
    float fuserate;             ///< rate of fusionreads in totalreads
    std::string gene1;          ///< hgene
    std::string chr1;           ///< chr on which hgene situated
    int32_t junctionposition1;  ///< junction position of hgene on chr1
    std::string strand1;        ///< strand hgene situated on chr1
    std::string transcript1;    ///< transcript of hgene
    std::string gene2;          ///< tgene
    std::string chr2;           ///< chr on which tgene situated
    int32_t junctionposition2;  ///< junction position of tgene on chr2
    std::string strand2;        ///< strand tgene situated on chr2
    std::string transcript2;    ///< transcript of tgene
    std::string fusionsequence; ///< fusion concensus sequence
    int32_t fseqbp;             ///< fusion sequence length of one partner in fusiongene
    std::string indb;           ///< this fusion event is in public database if true
    std::string svt;            ///< svtype, raw str representation
    int32_t svsize;             ///< sv size of this fusion event
    int32_t srcount;            ///< seed split reads supporting this fusion event
    int32_t dpcount;            ///< seed discordant pair of reads supporting this fusion event
    int32_t srrescued;          ///< all split reads supporting this fusion event
    int32_t dprescued;          ///< all siscordant pair of reads supporting this fusion event
    int32_t srrefcount;         ///< all reads supporting ref type
    int32_t dprefcount;         ///< all paired reads supporting ref type
    int32_t insbp;              ///< inserted sequence length around breakpoint
    std::string insseq;         ///< inserted sequence
    int32_t svid;               ///< sv id
    int32_t svint;              ///< svtype, precise integer representation
    TFUSION_FLAG fsmask;        ///< fusion mask
    std::string ts1name;        ///< transcript1 name(RNA only)
    int32_t ts1pos;             ///< transcript1 breakpoint(RNA only)
    std::string ts2name;        ///< transcript2 name(RNA only)
    int32_t ts2pos;             ///< transcript2 breakpoint(RNA only)

    /** construct an FusionRecord */
    FusionRecord(){
        fsmask = 0;
    }

    /** destroy an FusionRecord */
    ~FusionRecord(){}

    /** output an FusionRecord to ostream
     * @param os reference of ostream
     * @param fsr reference of FusionRecord
     * @return reference of ostream
     */
    inline friend std::ostream& operator<<(std::ostream& os, const FusionRecord& fsr){
        os << fsr.fusegene << "\t" << fsr.fusepattern << "\t";
        os << fsr.fusionreads << "\t" << fsr.totalreads << "\t" << fsr.fuserate << "\t";
        os << fsr.gene1 << "\t" << fsr.chr1 << "\t" << fsr.junctionposition1 << "\t" << fsr.strand1 << "\t" << fsr.transcript1 << "\t";
        os << fsr.gene2 << "\t" << fsr.chr2 << "\t" << fsr.junctionposition2 << "\t" << fsr.strand2 << "\t" << fsr.transcript2 << "\t";
        os << fsr.fusionsequence << "\t" << fsr.fseqbp << "\t" << fsr.indb << "\t" << fsr.svt << "\t" << fsr.svsize << "\t";
        os << fsr.srcount << "\t" << fsr.dpcount << "\t" << fsr.srrescued << "\t" << fsr.dprescued << "\t";
        os << fsr.srrefcount << "\t" << fsr.dprefcount << "\t";
        os << fsr.insbp << "\t" << fsr.insseq << "\t" << fsr.svid << "\t" << fsr.svint << "\t" << fsr.fsmask;
        if(fsr.fsmask & FUSION_FCALLFROMRNASEQ){
            os << "\t" << fsr.ts1name << "\t" << fsr.ts1pos << "\t" << fsr.ts2name << "\t" << fsr.ts2pos;
        }
        os << "\n";
       return os;
    }
};

/** class to store a list of FusionRecords */
typedef std::vector<FusionRecord> FusionRecordList;

#endif
