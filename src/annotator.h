#ifndef ANNOTATOR_H
#define ANNOTATOR_H

#include <util.h>
#include "stats.h"
#include "resbamstat.h"
#include <htslib/tbx.h>
#include <htslib/kstring.h>

/** SV breakpoint coverage annotator */
class Annotator{
    public:
        Options* mOpt; ///< pointer to Options

    public:
        /** default constructor of Annotator */
        Annotator(){}

        /** constructor of Annotator
         * @param opt pointer to Options object
         */
        Annotator(Options* opt) : mOpt(opt) {}
        
        /** destructor of Annotator */
        ~Annotator(){}

        /** refine coverage stat */
        void refineCovAnno(Stats* sts, const SVSet& svs);

        /** annotate SV coverage
         * @param svs reference of SVRecords
         */
        Stats* covAnnotate(SVSet& svs);

        /** split a cgranges_t by contig and label sum
         * @param cr pointer to cgranget_st
         * @param us unit sum to split
         * @param re result list
         */
        void cgrsplit(const cgranges_t* cr, std::vector<RegItemCnt>& ctgRng, int32_t us = 1000);

        /** annotate DNA SV gene information
         * @param svs reference of SVRecords
         * @param gl reference of GeneInfoList
         */
        void geneAnnoDNA(SVSet& svs, GeneInfoList& gl);

        /** annotate DNA SV gene information in a range of svs
         * @param svs reference of SVRecords
         * @param gl reference of GeneInfoList
         * @param begIdx beginning index of svs to annotate
         * @param endIdx ending index of svs to annotate
         */
        void rangeGeneAnnoDNA(SVSet& svs, GeneInfoList& gl, int32_t begIdx, int32_t endIdx);
        
        /** annotate RNA SV gene information
         * @param svs reference of SVRecords
         * @param gl reference of GeneInfoList
         */
        void geneAnnoRNA(SVSet& svs, GeneInfoList& gl);

        /** annotate RNA SV gene information in a range of svs
         * @param svs reference of SVRecords
         * @param gl reference of GeneInfoList
         * @param begIdx beginning index of svs to annotate
         * @param endIdx ending index of svs to annotate
         */
        void rangeGeneAnnoRNA(SVSet& svs, GeneInfoList& gl, int32_t begIdx, int32_t endIdx);

        /** get first overlap of an region against an region set
         * @param s region sets
         * @param p region to query
         * @return iterator to an region overlap with p or end
         */
        static std::set<std::pair<int32_t, int32_t>>::iterator getFirstOverlap(const std::set<std::pair<int32_t, int32_t>>& s, const std::pair<int32_t, int32_t>& p);

        /** get transcripts spanning a genome breakpoint
         * @param trl transcript list spanning this bp
         * @param chr chr of bp
         * @param pos pos of bp
         * @param fp file pointer of dna annotation db file
         * @param tbx index pointer of dna annotation db file
         */
        static void getDNABpTrs(TrsRecList& trl, const std::string& chr, int32_t pos, htsFile* fp, tbx_t* tbx);
        
        /** get transcripts spanning a transcriptnome breakpoint
         * @param trl transcript list spanning this bp
         * @param chr chr of bp
         * @param pos pos of bp
         * @param fp file pointer of dna annotation db file
         * @param tbx index pointer of dna annotation db file
         * @param isbp1 this is the starting position of sv if true
         * @param svt sv type of this sv event
         */
        static void getRNABpTrs(TrsRecList& trl, const std::string& chr, int32_t pos, htsFile* fp, tbx_t* tbx, bool isbp1, int32_t svt);
};

#endif
